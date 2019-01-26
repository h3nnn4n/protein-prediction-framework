import rosetta_pack
import protein_data
import random
import numpy as np
import math
import sys
import time
import datetime
import string

from operators.operators import Operators


class DE:
    def __init__(self, pop_size=50, pname='1zdd', c_rate=0.5, f_factor=0.5, max_iters=100):
        self.rosetta_pack = rosetta_pack.RosettaPack(pname)

        self.pname = pname

        self.pop_size = pop_size
        self.pop = [protein_data.ProteinData(self.rosetta_pack) for _ in range(pop_size)]
        self.c_rate = c_rate
        self.f_factor = f_factor

        self.trial = protein_data.ProteinData(self.rosetta_pack)
        self.max_iters = max_iters
        self.spent_iters = 0
        self.m_nmdf = 0
        self.cname = ''

        self.d = self.pop[0].pose.total_residue()

        # Energy function
        self.energy_function = None
        self.current_energy_function = None
        self.energy_options = []
        self.parsed_energy_options = {}
        self.spent_gens = 0
        self.spent_on_score = 0
        self.total_evals_on_current_score = 0

        # ETA
        self.buffer_size = 10
        self.time_buffer = [0 for _ in range(self.buffer_size)]
        self.spent_eval_buffer = [0 for _ in range(self.buffer_size)]
        self.last_time = 0
        self.last_spent_evals = 0
        self.time_pivot = 0
        self.start_time = 0
        self.end_time = 0
        self.repack_time = 0

        # Other stuff
        self.extended_diversity_measurements = False
        self.coil_only = False

        self.mode = None
        self.stop_condition = ''

        self.failsafe_verbose = False

        self.stage0_init = False
        self.stage0_init = False
        self.stage2_interval = -1
        self.stage2_all_interval = -1
        self.partial_reset = -1
        self.log_interval = 10

        self.reset_d_trigger = 0.0
        self.reset_d_percent = 0.0

        self.reset_rmsd_trigger = 0.0
        self.reset_rmsd_percent = 0.0

        self.config_name = None
        self.stats_name = None
        self.ops_stats_name = None
        self.stats = None
        self.ops_stats = None

        self.spent_evals = 0
        self.max_evals = 500000

        # Pop data dump
        # nothing

        # Clearing
        self.do_clearing = False
        self.clearing_interval = 10
        self.clearing_size = 1

        # REMC
        self.enable_remc = True
        self.update_remc()

        # Crowding
        self.do_crowding = False
        self.do_rmsd_crowding = False
        self.crowding_factor = None

        # Moment of Inertia
        self.centroids = None
        self.moment_of_inertia = None

        # LHS parameters
        self.do_lsh = False
        self.n_hashes = None
        self.n_buckets = None
        self.hashes = None
        self.hash_values = None
        self.active_hash1 = 0
        self.active_hash2 = 1
        self.tmp1 = None
        self.tmp2 = None
        self.update_interval = None
        self.change_interval = None

        # Log stuff
        self.init_time = time.time()
        self.now = datetime.datetime.now()

        # SaDE stuff
        self.sade_run = True
        self.sade_lp = 50
        self.sade_lp_left = self.sade_lp
        # self.sade_f = []
        self.sade_selection = ''

        self.ops = []

        self.sade_n_ops = None
        self.sade_ops_probs = None

        self.sade_success_memory = None
        self.sade_failure_memory = None

        self.sade_cr = None
        self.sade_cr_m = None
        self.sade_cr_memory = None

        self.sade_reinit_interval = None

        # Inner info
        self.mean = 0
        self.best_index = 0
        self.best_score = float('inf')

        self.mean_last_improv = None
        self.mean_improv_value = 0
        self.mean_last_improv_value = 0
        self.mean_improv_threshold = 0.1
        self.mean_improv_iter_threshold = 1000

        self.last_improv = None
        self.improv_value = 0
        self.last_improv_value = 0
        self.improv_threshold = 0.1
        self.improv_iter_threshold = 5000

        self.operators = Operators(de=self)

        print('Finished initialization')

# ######################### Config Update ##########################

    def reload_config(self):
        self.update_remc()

    def update_remc(self):
        self.trial.enable_remc = self.enable_remc

        for p in self.pop:
            p.enable_remc = self.enable_remc

# ######################### START OF SADE ##########################

    def set_sade_ops(self):
        self.sade_ops = []

        for op in self.ops:
            if op in self.operators.available_operators.keys():
                self.sade_ops.append(self.operators.get_operator(op))
            else:
                raise ValueError('An invalid operator was requested: %s' % op)

        self.sade_n_ops = len(self.sade_ops)

    def sade_reinit(self):
        print('reinit')
        self.sade_cr = [[random.random() for k in range(self.sade_n_ops)] for i in range(self.pop_size)]
        self.sade_cr_m = [np.clip(random.gauss(0.75, 0.1), 0.6, 1.0) for k in range(self.sade_n_ops)]
        self.sade_cr_memory = [[self.sade_cr_m[k]] for k in range(self.sade_n_ops)]

        self.sade_success_memory = [[0 for k in range(self.sade_n_ops)] for i in range(self.sade_lp)]
        self.sade_failure_memory = [[0 for k in range(self.sade_n_ops)] for i in range(self.sade_lp)]

        self.sade_ops_probs = [1 / self.sade_n_ops for _ in range(self.sade_n_ops)]

        self.sade_lp_left = self.sade_lp

    def sade_update_parameters(self):
        if self.sade_lp_left <= 0:
            for k in range(self.sade_n_ops):
                self.sade_cr_memory[k].sort()

            for k in range(self.sade_n_ops):
                if len(self.sade_cr_memory[k]) == 0:
                    print('Empty memory for k = %d %s!' % (k, self.sade_ops[k]))
                else:
                    self.sade_cr_m[k] = np.median(self.sade_cr_memory[k])

            self.sade_cr = [[np.clip(random.gauss(self.sade_cr_m[k], 0.1), 0.0, 1.0) for k in range(self.sade_n_ops)]
                            for i in range(self.pop_size)]
        else:
            self.sade_lp_left -= 1

    def sade_update_ops(self):
        if self.sade_lp_left <= 0:
            if self.sade_selection == 'roulette':
                for k in range(self.sade_n_ops):
                    s_s = sum([self.sade_success_memory[i][k] for i in range(self.sade_lp)])
                    s_f = sum([self.sade_failure_memory[i][k] for i in range(self.sade_lp)])
                    if s_s + s_f > 0:
                        self.sade_ops_probs[k] = s_s / (s_s + s_f) + 0.01
                    else:
                        self.sade_ops_probs[k] = 0.01
            elif 'tournament' in self.sade_selection:
                for k in range(self.sade_n_ops):
                    s_s = sum([self.sade_success_memory[i][k] for i in range(self.sade_lp)])
                    s_f = sum([self.sade_failure_memory[i][k] for i in range(self.sade_lp)])
                    if s_s - s_f > 0:
                        self.sade_ops_probs[k] = (s_s - s_f) + 0.01
                    else:
                        self.sade_ops_probs[k] = 0.01

        norm = sum(self.sade_ops_probs)
        self.sade_ops_probs = list(map(lambda x: x / norm, self.sade_ops_probs))

    def sade_get_op(self):
        if self.sade_selection == 'roulette':
            n = random.random()
            a = 0.0
            i = -1

            while a < n:
                i += 1
                a += self.sade_ops_probs[i]
        elif 'tournament' in self.sade_selection:
            if self.sade_selection == 'tournament':
                size = 2
            else:
                size = int(self.sade_selection[10:])

            trial = random.sample(list(range(self.sade_n_ops)), k=size)
            t_values = [self.sade_ops_probs[i] for i in trial]

            w = 0
            i = None
            for k, v in zip(trial, t_values):
                if v > w:
                    w = v
                    i = k

            if i is None:
                print('rip')
                print(trial, t_values)
                print(self.sade_ops_probs)

        return self.sade_ops[i]

# ######################### END OF SADE ##########################

    def open_stats(self):
        char_set = string.ascii_uppercase + string.digits
        r_string = ''.join(random.sample(char_set * 6, 6))

        now = self.now

        self.name_suffix = "_%s__%s__%04d_%02d_%02d__%02d_%02d_%02d__%s" % \
                           (self.pname, self.cname, now.year, now.month, now.day, now.hour, now.minute, now.second, r_string)
        self.stats_name = self.rosetta_pack.protein_loader.original + '/' + "stats_" + self.name_suffix + ".dat"
        self.stats = open(self.stats_name, 'w')

        self.ops_stats_name = self.rosetta_pack.protein_loader.original + '/' + "ops_" + self.name_suffix + ".dat"
        self.ops_stats = open(self.ops_stats_name, 'w')

    def dump_config(self):
        self.config_name = self.rosetta_pack.protein_loader.original + '/' + "parameters_" + self.name_suffix + ".yaml"
        with open(self.config_name, 'w') as f:
            f.write('pname: %s\n' % (self.pname))
            f.write('pop_size: %d\n' % (self.pop_size))
            f.write('c_rate: %f\n' % (self.c_rate))
            f.write('f_factor: %f\n' % (self.f_factor))
            f.write('max_iters: %d\n' % (self.max_iters))
            f.write('max_evals: %d\n' % (self.max_evals))
            f.write('stop_condition: %s\n' % (self.stop_condition))
            f.write('stage0_init: %d\n' % self.stage0_init)
            f.write('stage2_interval: %d\n' % self.stage2_interval)
            f.write('stage2_all_interval: %d\n' % self.stage2_all_interval)
            f.write('partial_reset: %d\n' % self.partial_reset)
            f.write('log_interval: %d\n' % self.log_interval)
            f.write('island_interval: %d\n' % self.island_interval)
            f.write('do_lsh: %d\n' % self.do_lsh)
            f.write('n_hashes: %d\n' % self.n_hashes)
            f.write('update_interval: %d\n' % self.update_interval)
            f.write('change_interval: %d\n' % self.change_interval)
            f.write('reset_d_trigger: %f\n' % self.reset_d_trigger)
            f.write('reset_d_percent: %f\n' % self.reset_d_percent)
            f.write('reset_rmsd_percent: %f\n' % self.reset_rmsd_percent)
            f.write('reset_rmsd_trigger: %f\n' % self.reset_rmsd_trigger)
            f.write('sade_run: %d\n' % self.sade_run)
            f.write('sade_lp: %d\n' % self.sade_lp)
            f.write('sade_reinit_interval: %d\n' % self. sade_reinit_interval)
            f.write('sade_selection: %s\n' % self.sade_selection)
            f.write('enable_remc: %s\n' % self.enable_remc)
            f.write('do_crowding: %d\n' % self.do_crowding)
            f.write('do_rmsd_crowding: %d\n' % self.do_rmsd_crowding)
            f.write('crowding_factor: %d\n' % self.crowding_factor)
            f.write('ops: %s\n' % self.ops)
            f.write('do_clearing: %d\n' % self.do_clearing)
            f.write('clearing_interval: %d\n' % self.clearing_interval)
            f.write('clearing_size: %d\n' % self.clearing_size)
            f.write('energy_function: %s\n' % self.energy_function)
            f.write('energy_options: %s\n' % self.energy_options)
            f.write('extended_diversity_measurements: %s\n' % self.extended_diversity_measurements)
            f.flush()

    def set_allatom(self):
        switch = self.rosetta_pack.allatom_switch
        pack = self.rosetta_pack.get_packer()
        mini = self.rosetta_pack.get_min_mover()

        for p in self.pop:
            switch.apply(p.pose)
            pack.apply(p.pose)
            mini.apply(p.pose)
            p.update_angle_from_pose()
            p.eval()

# ######################### START OF LSH ##########################

    def create_hashs(self):
        self.hashes = [np.random.randint(100, size=self.pop[0].nsca) for _ in range(self.n_hashes)]

    def apply_hash(self, debug=False):
        # debug = True

        if self.hash_values is None:
            self.hash_values = [[] for _ in range((self.n_buckets + 1) ** 2 + (self.n_buckets + 1))]
        else:
            for k, v in enumerate(self.hash_values):
                if len(v) > 0:
                    self.hash_values[k] = []

        minh1 = None
        maxh1 = None
        minh2 = None
        maxh2 = None

        if self.tmp1 is None:
            self.tmp1 = np.empty(self.pop_size)

        if self.tmp2 is None:
            self.tmp2 = np.empty(self.pop_size)

        tmp1 = self.tmp1
        tmp2 = self.tmp2

        for i in range(self.pop_size):
            h1 = np.dot(self.hashes[self.active_hash1], self.pop[i].angles)
            h2 = np.dot(self.hashes[self.active_hash2], self.pop[i].angles)
            tmp1[i] = h1
            tmp2[i] = h2

            if minh1 is None or h1 < minh1:
                minh1 = h1

            if maxh1 is None or h1 > maxh1:
                maxh1 = h1

            if minh2 is None or h2 < minh2:
                minh2 = h2

            if maxh2 is None or h2 > maxh2:
                maxh2 = h2

        r1 = (maxh1 - minh1) / self.n_buckets
        b1 = random.random() * r1

        r2 = (maxh2 - minh2) / self.n_buckets
        b2 = random.random() * r2

        if debug:
            print("r1: %8.3f  b1: %8.3f  minh1: %8.3f  maxh1: %8.3f" % (r1, b1, minh1, maxh1))
            print("r2: %8.3f  b2: %8.3f  minh2: %8.3f  maxh2: %8.3f" % (r2, b2, minh2, maxh2))

        for i in range(self.pop_size):
            a, b = math.floor((tmp1[i] - minh1 + b1) / r1), math.floor((tmp2[i] - minh2 + b2) / r2)
            v = (self.n_buckets) * a + b

            try:
                self.hash_values[v].append(i)
            except Exception:
                raise IndexError(v, len(self.hash_values))

        for n, i in enumerate(self.hash_values):
            if debug:
                if len(i) > 0:
                    print(n, i)

        if debug:
            print()

    def change_hash(self):
        old1 = self.active_hash1
        old2 = self.active_hash2

        while old1 == self.active_hash1:
            self.active_hash1 = random.randint(0, self.n_hashes - 1)

        while old2 == self.active_hash2 and self.active_hash1 == self.active_hash2:
            self.active_hash2 = random.randint(0, self.n_hashes - 1)

# ######################### END OF LSH ##########################

    def get_best(self):
        return self.pop[self.best_index].angles

    def run(self):
        self.set_sade_ops()
        self.open_stats()
        self.dump_config()

        self.create_hashs()

        self.update_threshold()
        self.update_score_function(step=False)
        # self.update_mean()

        if self.stage0_init:
            print('Stage0 init')
            # self.log(it=-1)
            for p in self.pop:
                p.eval()
                p.stage1_mc()
                p.update_angle_from_pose()
                p.eval()
            self.trial.eval()
        else:
            for p in self.pop:
                p.eval()
            self.trial.eval()

        self.update_mean()

        # self.log(it=0)

        if self.do_lsh:
            self.apply_hash()

        self.sade_reinit()
        self.start_time = time.time()
        self.last_time = self.start_time

        self.stop_at_eval = 'evals' in self.stop_condition
        self.stop_at_iter = 'iters' in self.stop_condition

        self.it = 0
        while (self.spent_evals < self.max_evals and self.stop_at_eval) or (self.it < self.max_iters and self.stop_at_iter) or self.mode == 'marathon':
            if self.sade_run and (self.sade_reinit_interval > 0 and self.it % self.sade_reinit_interval == 0) and \
               self.energy_function not in ['cascade']:
                self.sade_reinit()

            if self.sade_run:
                self.sade_update_parameters()
                self.sade_update_ops()

            self.it += 1

            self.update_score_function()

            if self.do_lsh and self.it % self.change_interval == 0:
                self.change_hash()

            if self.do_lsh and self.it % self.update_interval == 0:
                self.apply_hash()

            for i in range(self.pop_size):
                self.huehue = i
                if self.do_lsh and not self.sade_run:
                    ret = self.rand1bin_lsh(i)  # pylint: disable=E1128
                else:
                    if self.sade_run:
                        ret = self.sade_get_op()(i)
                    else:
                        ret = self.operators.rand1bin_global(i)  # pylint: disable=E1111

                if ret is None:
                    self.spent_evals += 1
                else:
                    self.spent_evals += ret

            if self.do_lsh and False:
                for h in self.hash_values:
                    if len(h) >= self.pop_size // 2:
                        print('Niche reset')
                        for i in h:
                            if random.random() < .5 and i != self.best_index:
                                self.pop[i].stage1_mc()
                                self.pop[i].update_angle_from_pose()
                                self.pop[i].eval()
                                self.pop[i].stage2_mc()
                                self.pop[i].update_angle_from_pose()
                                self.pop[i].eval()

            if self.it % 50 == 0 and self.avg_rmsd() < self.reset_rmsd_trigger:
                print('rmsd_reset')
                # self.log()
                # self.trial.stage1_mc(n=100)
                # self.trial.reset()
                for i in range(self.pop_size):
                    if random.random() < self.reset_rmsd_percent and i != self.best_index:
                        # before = self.pop[i].score
                        # import ipdb; ipdb.set_trace()
                        self.pop[i].reset()
                        self.pop[i].stage1_mc(n=100)
                        # self.pop[i].update_angle_from_pose()
                        self.pop[i].eval()
                        # print(i, before, self.pop[i].score)
                self.update_diversity()
                # self.log()

            if self.it > 1 and self.diversity < self.reset_d_trigger:
                print('d_reset')
                for i in range(self.pop_size):
                    if random.random() < self.reset_d_percent and i != self.best_index:
                        self.pop[i].reset()
                        self.pop[i].stage1_mc()
                        self.pop[i].update_angle_from_pose()
                        self.pop[i].eval()
                self.update_diversity()

            if (self.partial_reset > 0 and self.it % self.partial_reset == 0 and self.it > 0):
                print('Partial reset')
                for i in range(self.pop_size):
                    if random.random() < .15 and i != self.best_index:
                        self.pop[i].reset()
                        self.pop[i].stage2_mc()
                        self.pop[i].update_angle_from_pose()
                        self.pop[i].eval()
                self.update_diversity()

            if self.stage2_interval > 0 and self.it % self.stage2_interval == 0 and self.it > 0:
                print('LS')
                self.pop[self.best_index].stage2_mc()
                self.pop[self.best_index].update_angle_from_pose()
                self.pop[self.best_index].eval()

            if self.stage2_all_interval > 0 and self.it % self.stage2_all_interval == 0 and self.it > 0:
                print('NINJA MOVE')
                self.rosetta_pack.loop_modeling(self.pop[self.best_index].pose)
                self.pop[self.best_index].update_angle_from_pose()
                self.pop[self.best_index].eval()
                for i in range(self.pop_size):
                    if i == self.best_index:
                        continue

                    if random.random() < .25:
                        self.pop[i].stage2_mc(n=10, temp=2.5)
                    else:
                        self.pop[i].stage2_mc(n=10, temp=.5)
                    self.pop[i].update_angle_from_pose()
                    self.pop[i].eval()

            self.update_mean()
            self.update_threshold()

            if False and self.it % 1000 == 0:
                self.dump_pbd_best(self.it)

            if self.log_interval > 0 and self.it % self.log_interval == 0 or self.it == 1:
                self.log()
                self.stats.flush()

                sys.stdout.flush()

            self.rosetta_pack.pymover.apply(self.pop[self.best_index].pose)

        self.end_time = time.time()

        self.log()
        self.dump_pbd_best(self.it)

        rmsd = self.rosetta_pack.get_rmsd_from_pose(self.pop[self.best_index].pose)
        oldscore = self.best_score
        score = self.pop[self.best_index].repack()

        name = self.rosetta_pack.protein_loader.original + '/' + ("best_repacked_%05d_" % self.it) + self.name_suffix + ".pdb"
        self.pop[self.best_index].repacked.dump_pdb(name)

        repack_name = name

        tm_before = self.pop[self.best_index].run_tmscore()
        self.rosetta_pack.run_tmscore(name=repack_name)
        tm_after = self.rosetta_pack.get_tmscore()

        self.repack_time = time.time()
        name = self.rosetta_pack.protein_loader.original + '/' + "repack_" + self.name_suffix + ".dat"
        with open(name, 'w') as f:
            f.write('repack_time:        %12.4f\n' % (self.repack_time - self.end_time))
            f.write('score:              %12.4f\n' % oldscore)
            f.write('scorefxn:           %12.4f\n' % score)
            f.write('rmsd_after:         %12.4f\n' % (self.rosetta_pack.get_rmsd_from_pose(self.pop[self.best_index].repacked)))
            f.write('rmsd_before:        %12.4f\n' % rmsd)
            f.write('rmsd_change:        %12.4f\n' % (rmsd - self.rosetta_pack.get_rmsd_from_pose(self.pop[self.best_index].repacked)))
            f.write('tm_score_before:    %12.4f\n' % tm_before['tm_score'])
            f.write('maxsub_before:      %12.4f\n' % tm_before['maxsub'])
            f.write('gdt_ts_before:      %12.4f\n' % tm_before['gdt_ts'][0])
            f.write('gdt_ha_before:      %12.4f\n' % tm_before['gdt_ha'][0])
            f.write('tm_score_after:     %12.4f\n' % tm_after['tm_score'])
            f.write('maxsub_after:       %12.4f\n' % tm_after['maxsub'])
            f.write('gdt_ts_after:       %12.4f\n' % tm_after['gdt_ts'][0])
            f.write('gdt_ha_after:       %12.4f\n' % tm_after['gdt_ha'][0])
            f.write('tm_score_change:    %12.4f\n' % (tm_before['tm_score'] - tm_after['tm_score']))
            f.write('maxsub_change:      %12.4f\n' % (tm_before['maxsub'] - tm_after['maxsub']))
            f.write('gdt_ts_change:      %12.4f\n' % (tm_before['gdt_ts'][0] - tm_after['gdt_ts'][0]))
            f.write('gdt_ha_change:      %12.4f\n' % (tm_before['gdt_ha'][0] - tm_after['gdt_ha'][0]))
            f.write('gdt_ts_info_before: %12.4f %12.4f %12.4f %12.4f\n' % (tm_before['gdt_ts'][1][0], tm_before['gdt_ts'][1][0], tm_before['gdt_ts'][1][0], tm_before['gdt_ts'][1][0]))
            f.write('gdt_ha_info_after:  %12.4f %12.4f %12.4f %12.4f\n' % (tm_after['gdt_ha'][1][0], tm_after['gdt_ha'][1][0], tm_after['gdt_ha'][1][0], tm_after['gdt_ha'][1][0]))

    def print_hash(self):
        for n, i in enumerate(self.hash_values):
            if len(i) > 0:
                print(n, i)

# ########### End of operators

    def selection(self, candidate=None):
        if self.do_crowding or self.do_rmsd_crowding:
            self.crowding_selection()
        else:
            if candidate is not None:
                self.standard_selection(candidate)
            else:
                pass

    def crowding_selection(self):
        sade_k = self.sade_k

        _, cr = self.get_f_cr()

        ps = random.sample(range(self.pop_size), k=self.crowding_factor)

        best_index = None
        best_dist = None

        for p in ps:
            if self.do_rmsd_crowding:
                d = self.rosetta_pack.get_rmsd_from_pose(self.trial.pose, self.pop[p].pose)
            elif self.do_crowding:
                d = np.sqrt(np.sum((self.trial.angles - self.pop[p].angles)**2))
            else:
                print('WARNING! Crowding was called without being enabled!')

            if best_dist is None or d < best_dist:
                best_dist = d
                best_index = p

        candidate = self.pop[best_index]

        if self.trial.score < candidate.score:
            if self.sade_run:
                self.sade_cr_memory[sade_k].append(cr)
                ind = self.it % self.sade_lp
                self.sade_success_memory[ind][sade_k] += 1

            t = candidate
            candidate = self.trial
            self.trial = t
            if self.trial is candidate:
                import sys
                print('PANIC! Found duplicated reference in population')
                sys.exit()
        else:
            if self.sade_run:
                ind = self.it % self.sade_lp
                self.sade_failure_memory[ind][sade_k] += 1

    def standard_selection(self, candidate):
        sade_k = self.sade_k

        _, cr = self.get_f_cr()

        if self.trial.score < candidate.score:
            if self.sade_run:
                self.sade_cr_memory[sade_k].append(cr)
                ind = self.it % self.sade_lp
                self.sade_success_memory[ind][sade_k] += 1

            for k in range(self.pop_size):
                p = self.pop[k]
                if candidate is p:
                    t = self.pop[k]
                    self.pop[k] = self.trial
                    self.trial = t
        else:
            if self.sade_run:
                ind = self.it % self.sade_lp
                self.sade_failure_memory[ind][sade_k] += 1

    def get_f_cr(self):
        if self.sade_run:
            return self.get_f(), self.get_cr()

        return self.f_factor, self.c_rate

    def get_f(self):
        return random.gauss(0.5, 0.3)

    def get_cr(self):
        return self.sade_cr[self.huehue][self.sade_k]

    def update_score_function(self, step=True):
        def update_pop_score(update_score=False):
            for p in self.pop:
                p.set_score_function(self.current_energy_function)
                if update_score:
                    p.eval()

            self.trial.set_score_function(self.current_energy_function)
            if update_score:
                self.trial.eval()
                self.update_mean()

        if self.energy_function in ['score0', 'score1', 'score2', 'score3', 'score5']:
            if not hasattr(self, 'single_score_update'):
                self.single_score_update = True
                self.current_energy_function = self.energy_function
                update_pop_score()

        elif self.energy_function == 'cascade':
            if self.spent_gens == 0 and not hasattr(self, 'parsed_energy'):
                first = self.parse_energy_options()
                self.parsed_energy = True
                self.current_energy_function = first
                update_pop_score()
            if self.mode == 'marathon':
                low, high, n = self.parsed_energy_options[self.current_energy_function]
                # print(self.spent_gens, high, low, self.current_energy_function, n)
                # if (self.best_score is not None and abs(self.mean - self.best_score) < 0.1):
                if self.last_improv > self.improv_iter_threshold or self.mean_last_improv > self.mean_improv_iter_threshold:
                    if self.last_improv > self.improv_iter_threshold:
                        print('Iter reset!')
                        # print('Next: last improv', self.last_improv, self.improv_iter_threshold)
                    elif self.mean_last_improv > self.mean_improv_iter_threshold:
                        print('Mean reset!')
                        # print('Next: last mean improv', self.mean_last_improv, self.mean_improv_iter_threshold)
                    if self.current_energy_function == n:
                        # import sys
                        # sys.exit()
                        self.mode = 'quit'
                    self.current_energy_function = n
                    update_pop_score(update_score=True)
                    self.sade_reinit()
                    self.spent_on_score = 0
                    self.total_evals_on_current_score = self.spent_on_score

                # print(self.spent_gens, self.spent_on_score, self.total_evals_on_current_score)
                if self.sade_reinit_interval > 0 and self.spent_on_score % self.sade_reinit_interval == 0 and self.sade_run and \
                   self.spent_on_score > 0 and self.spent_on_score != self.total_evals_on_current_score:
                    self.sade_reinit()

                if step:
                    self.spent_gens += 1
                    self.spent_on_score += 1
            else:
                low, high, n = self.parsed_energy_options[self.current_energy_function]
                # print(self.spent_gens, high, low, self.current_energy_function, n)
                if high < self.spent_gens:
                    self.current_energy_function = n
                    update_pop_score(update_score=True)
                    self.sade_reinit()
                    self.spent_on_score = high - low
                    self.total_evals_on_current_score = self.spent_on_score

                if self.spent_on_score % self.sade_reinit_interval == 0 and self.sade_run and self.spent_on_score > 0 and \
                   self.spent_on_score < self.total_evals_on_current_score:
                    self.sade_reinit()

                if step:
                    self.spent_gens += 1
                    self.spent_on_score -= 1

    def parse_energy_options(self):
        first = None
        acc = 0
        for w in self.energy_options:
            a, b = w.split('_')
            # if acc == 0:
            #     self.current_energy_function = a
            b = int(b)
            acc += b
            # print(a, b, acc - b, acc)
            self.parsed_energy_options[a] = (acc - b, acc)

            if first is None:
                first = a

        last = None
        for k, w in enumerate(self.energy_options):
            a, _ = w.split('_')

            if k == 0:
                last = a
                continue

            x, y = self.parsed_energy_options[last]
            # print(last, x, y)
            self.parsed_energy_options[last] = (x, y, a)
            # print(last, self.parsed_energy_options[last])
            last = a

        x, y = self.parsed_energy_options[last]
        self.parsed_energy_options[last] = (x, y, a)
        # print(last, self.parsed_energy_options[last])
        # print(last)

        return first

    def update_mean(self):
        self.mean = 0
        self.best_score = float('inf')
        for i in range(self.pop_size):
            self.mean += self.pop[i].score / self.pop_size
            if self.best_score is None or self.pop[i].score < self.best_score:
                self.best_score = self.pop[i].score
                self.best_index = i

        return self.mean

    def update_threshold(self, force_update=False):
        if self.mean_last_improv is None:
            self.mean_last_improv = 0
            self.mean_last_improv_value = float("inf")
        # elif abs(self.mean_last_improv_value - self.mean) > self.mean_improv_threshold or self.mean_last_improv > self.mean_improv_iter_threshold:
        elif abs(self.mean_last_improv_value - self.mean) > self.mean_improv_threshold:
            # print('Mean reset')
            self.mean_improv_value = self.mean - self.mean_last_improv_value
            self.mean_last_improv_value = self.mean
            self.mean_last_improv = 0
        else:
            self.mean_last_improv += 1

        if self.last_improv is None:
            self.last_improv = 0
            self.last_improv_value = float("inf")
        # elif abs(self.last_improv_value - self.best_score) > self.improv_threshold or self.last_improv > self.improv_iter_threshold:
        elif abs(self.last_improv_value - self.best_score) > self.improv_threshold:
            # print('Iter reset')
            self.improv_value = self.best_score - self.last_improv_value
            self.last_improv_value = self.best_score
            self.last_improv = 0
        else:
            self.last_improv += 1

    def update_moment_of_inertia(self):
        pop_size = self.pop_size
        n_dim = self.pop[0].nsca

        if self.centroids is None:
            self.centroids = [0 for _ in range(n_dim)]

        for i in range(n_dim):
            self.centroids[i] = 0.0
            for j in range(pop_size):
                self.centroids[i] += self.pop[j].angles[i] / pop_size

        self.moment_of_inertia = 0.0

        for i in range(n_dim):
            w = 0
            for j in range(pop_size):
                w += (self.pop[j].angles[i] - self.centroids[i]) ** 2.0

            self.moment_of_inertia += math.sqrt(w) / (pop_size - 1)

        self.moment_of_inertia /= n_dim

        return self.moment_of_inertia

    def update_diversity(self):
        diversity = 0
        aux_1 = 0
        aux_2 = 0

        d = 0

        nsca = len(self.pop[0].angles)

        for i in range(0, self.pop_size):
            for j in range(i + 1, self.pop_size):
                aux_1 = 0

                ind_a = self.pop[i].angles
                ind_b = self.pop[j].angles

                for d in range(0, nsca):
                    aux_1 += (ind_a[d] - ind_b[d]) ** 2

                aux_1 = math.sqrt(aux_1) / nsca

                if j == i + 1 or aux_2 > aux_1:
                    aux_2 = aux_1

            diversity += math.log(1.0 + aux_2)

        self.m_nmdf = max(self.m_nmdf, diversity)

        if self.m_nmdf > 0:
            self.diversity = diversity / self.m_nmdf
        else:
            self.diversity = 0.0

        return self.diversity

    def dump_pop_data(self):
        data = []
        for p in self.pop:
            r = self.rosetta_pack.get_rmsd_from_pose(p.pose)
            s = p.score
            data.append((r, s))

        name = self.rosetta_pack.protein_loader.original + '/' + ("popdata_%08d_" % self.it) + self.name_suffix + ".dat"
        with open(name, 'w') as f:
            for a, b in data:
                f.write('%12.5f %12.5f\n' % (a, b))

    def dump_pbd_pop(self):
        pass

    def dump_pbd_best(self, it):
        name = self.rosetta_pack.protein_loader.original + '/' + ("best_%05d_" % self.it) + self.name_suffix + ".pdb"
        self.pop[self.best_index].pose.dump_pdb(name)

    def avg_distance(self):
        if not self.extended_diversity_measurements:
            return 0

        s = 0
        c = 0
        pop = self.pop

        for i in range(self.pop_size):
            for j in range(i + 1, self.pop_size):
                c += 1
                diff = pop[i].angles - pop[j].angles
                s += np.linalg.norm(diff)

        return s / c

    def avg_rmsd(self):
        if not self.extended_diversity_measurements:
            return 0

        s = 0
        c = 0

        for i in range(self.pop_size):
            for j in range(i + 1, self.pop_size):
                c += 1
                s += self.rosetta_pack.get_rmsd_from_pose(self.pop[i].pose, self.pop[j].pose)

        return s / c

    def avg_rmsd_from_native(self):
        if not self.extended_diversity_measurements:
            return 0

        s = 0
        c = 0

        for i in range(self.pop_size):
            c += 1
            s += self.rosetta_pack.get_rmsd_from_pose(self.pop[i].pose)

        return s / c

    def log(self, it=None):
        if it is None:
            it = self.it

        rmsd = self.rosetta_pack.get_rmsd_from_pose(self.pop[self.best_index].pose)

        if it < 1:
            secs_per_iter = 0
            secs_per_eval = 0
            eta_evals = 0
            eta_iters = 0
        else:
            now = time.time()
            dt = now - self.last_time
            de = (self.spent_evals - self.last_spent_evals)

            self.spent_eval_buffer[self.time_pivot % self.buffer_size] = de

            if self.time_pivot >= self.buffer_size:
                self.time_buffer[self.time_pivot % self.buffer_size] = dt
            else:
                for i in range(self.time_pivot, self.buffer_size):
                    self.time_buffer[i] = dt

            secs_per_iter = np.mean(self.time_buffer) / self.log_interval
            secs_per_eval = sum(self.time_buffer) / sum(self.spent_eval_buffer)
            eta_evals = (self.max_evals - self.spent_evals) * secs_per_eval

            if self.mode == 'marathon':
                eta_iters = self.max_iters * secs_per_iter
            else:
                eta_iters = (self.max_iters - it) * secs_per_iter

            self.last_time = now
            self.last_spent_evals = self.spent_evals
            self.time_pivot += 1

        cr = ''  # '%3.2f' % self.c_rate
        probs = ''

        if self.sade_run:
            if self.sade_cr_m is not None:
                for i in range(len(self.sade_cr_m)):
                    cr += '%3.2f ' % self.sade_cr_m[i]

            if self.sade_ops_probs is not None:
                for k in range(self.sade_n_ops):
                    probs += '%3.2f ' % self.sade_ops_probs[k]

        string = ''

        data = [
            ('%8d', self.spent_evals),
            ('%8d', it),
            ('%8.4f', self.best_score),
            ('%8.4f', self.mean),
            ('%8.4f', self.update_diversity()),
            ('%8.4f', self.avg_rmsd()),
            ('%8.4f', rmsd),
            ('%8.4f', self.update_moment_of_inertia()),
            ('%10.2f', (time.time() - self.start_time)),
            ('%8.5f', secs_per_eval),
            ('%8.2f', eta_evals),
            ('%8.5f', secs_per_iter),
            ('%8.2f', eta_iters),
            ('%s', cr),
            ('%s', probs)
        ]

        # print(self.sade_success_memory)
        # print(self.sade_failure_memory)

        for k, (a, b) in enumerate(data):
            try:
                string += a % b
                string += ' '
            except Exception:
                print(k, a, b)

        print(string)
        self.stats.write(string + '\n')

        # print("%8.4f %8d %8.4f %8d" % (self.improv_value,
        #                                self.improv_iter_threshold - self.last_improv,
        #                                self.mean_improv_value,
        #                                self.mean_improv_iter_threshold - self.mean_last_improv))

        if self.ops_stats is not None and self.sade_success_memory is not None:
            string = ''

            string += '%8d  ' % it

            for k in range(self.sade_n_ops):
                s_s = sum([self.sade_success_memory[i][k] for i in range(self.sade_lp)])
                s_f = sum([self.sade_failure_memory[i][k] for i in range(self.sade_lp)])
                string += '  %8d %8d %8d' % (s_s + s_f, s_s, s_f)
            self.ops_stats.write(string + '\n')
            self.ops_stats.flush()
