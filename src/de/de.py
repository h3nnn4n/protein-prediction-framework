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
from locality_sensitive_hashing import LocalitySensitiveHashing
from forced_insertion import ForcedInsertion
from piecewise_exchange import PiecewiseExchange
from hooke_jeeves_postprocessing import HookeJeevesPostprocessing
from repacker import Repacker


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

        self.stop_condition = ''

        self.failsafe_verbose = False

        self.stage0_init = False
        self.log_interval = 10

        # Config
        self.config_name = None
        self.stats_name = None
        self.ops_stats_name = None
        self.stats = None
        self.ops_stats = None

        self.spent_evals = 0
        self.max_evals = 500000

        # Forced insertion
        self.forced_insertion = False
        self.forced_insertion_chance = 0.0
        self.forced_insertion_mode = None

        # Piecewise Exchange
        self.piecewise_exchange = PiecewiseExchange(de=self)

        # REMC
        self.enable_remc = True
        self.update_remc()

        # LS
        self.run_hooke_jeeves_postprocessing = True
        self.hooke_jeeves_postprocessing_mode = 'best'
        self.hooke_jeeves_postprocessing = HookeJeevesPostprocessing(de=self)

        # Moment of Inertia
        self.centroids = None
        self.moment_of_inertia = None

        # LHS parameters
        self.do_lsh = False
        self.n_hashes = None
        self.n_buckets = None
        self.update_interval = None
        self.change_interval = None

        # Log stuff
        self.init_time = time.time()
        self.now = datetime.datetime.now()

        # SaDE stuff
        self.sade_run = True
        self.sade_lp = 50
        self.sade_lp_left = self.sade_lp
        self.sade_selection = ''
        self.sade_k = None

        self.ops = []

        self.sade_n_ops = None
        self.sade_ops_probs = None

        self.sade_success_memory = None
        self.sade_failure_memory = None

        self.sade_cr = None
        self.sade_cr_m = None
        self.sade_cr_memory = None

        self.sade_reinit_interval = None

        # Repacking
        self.repack_mode = 'best'
        self.repacker = Repacker(de=self)

        # Inner info
        self.mean = 0
        self.best_index = 0
        self.best_score = float('inf')

        self.operators = Operators(de=self)

        print('Finished initialization')

# ######################### Config Update ##########################

    def reload_config(self):
        self.update_remc()
        self.update_lsh()
        self.forced_insertion_op = ForcedInsertion(
            de=self,
            mode=self.forced_insertion_mode,
            chance=self.forced_insertion_chance,
            enable=self.forced_insertion
        )

    def update_remc(self):
        self.trial.enable_remc = self.enable_remc

        for p in self.pop:
            p.enable_remc = self.enable_remc

    def update_lsh(self):
        if not self.do_lsh:
            return

        self.locality_sensitive_hashing = LocalitySensitiveHashing(de=self)
        self.locality_sensitive_hashing.inject_parameters()

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

    def check_sade_reinit(self):
        reached_reinit_step = (
            self.sade_reinit_interval > 0 and
            self.it % self.sade_reinit_interval == 0
        )

        if self.sade_run and reached_reinit_step:
            self.sade_reinit()

    def sade_update(self):
        if self.sade_run:
            self.sade_update_parameters()
            self.sade_update_ops()

# ######################### Stage0 Init ##########################

    def stage0_pop_init(self):
        if self.stage0_init:
            print('Stage0 init')
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

# ######################### LOGGING ##############################

    def generation_logger(self):
        if self.log_interval > 0 and self.it % self.log_interval == 0 or self.it == 1:
            self.log()
            self.stats.flush()
            sys.stdout.flush()

# ######################### RUN ##################################

    def open_stats(self):
        char_set = string.ascii_uppercase + string.digits
        r_string = ''.join(random.sample(char_set * 6, 6))

        now = self.now

        self.name_suffix = "_%s__%s__%04d_%02d_%02d__%02d_%02d_%02d__%s" % (
            self.pname, self.cname, now.year, now.month, now.day, now.hour, now.minute, now.second, r_string)
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
            f.write('log_interval: %d\n' % self.log_interval)
            f.write('do_lsh: %d\n' % self.do_lsh)
            f.write('n_hashes: %d\n' % self.n_hashes)
            f.write('n_buckets: %d\n' % self.n_buckets)
            f.write('update_interval: %d\n' % self.update_interval)
            f.write('change_interval: %d\n' % self.change_interval)
            f.write('forced_insertion: %d\n' % self.forced_insertion)
            f.write('forced_insertion_chance: %f\n' % self.forced_insertion_chance)
            f.write('forced_insertion_mode: %s\n' % self.forced_insertion_mode)
            f.write('sade_run: %d\n' % self.sade_run)
            f.write('sade_lp: %d\n' % self.sade_lp)
            f.write('sade_reinit_interval: %d\n' % self. sade_reinit_interval)
            f.write('sade_selection: %s\n' % self.sade_selection)
            f.write('enable_remc: %s\n' % self.enable_remc)
            f.write('run_hooke_jeeves_postprocessing: %s\n' % self.run_hooke_jeeves_postprocessing)
            f.write('hooke_jeeves_postprocessing_mode: %s\n' % self.hooke_jeeves_postprocessing_mode)
            f.write('ops: %s\n' % self.ops)
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

    def get_best(self):
        return self.pop[self.best_index].angles

    def run(self):
        self.set_sade_ops()
        self.open_stats()
        self.dump_config()
        self.forced_insertion_op.initialize_logger()

        if self.do_lsh:
            self.locality_sensitive_hashing.create_hashs()

        self.stage0_pop_init()

        self.update_mean()

        if self.do_lsh:
            self.locality_sensitive_hashing.apply_hash()

        self.sade_reinit()
        self.start_time = time.time()
        self.last_time = self.start_time

        self.stop_at_eval = 'evals' in self.stop_condition
        self.stop_at_iter = 'iters' in self.stop_condition

        self.it = 0

        while self.can_run_evals() or self.can_run_generations():
            self.it += 1

            self.check_sade_reinit()
            self.sade_update()
            self.forced_insertion_op.run()

            if self.do_lsh:
                self.locality_sensitive_hashing.lsh_step()

            self.population_update()
            self.update_mean()
            self.generation_logger()
            self.rosetta_pack.pymover.apply(self.pop[self.best_index].pose)

        self.end_time = time.time()

        self.hooke_jeeves_postprocessing.run_hooke_jeeves()

        self.log()
        self.dump_pbd_best(self.it)

        self.repacker.run_repack()

# ######################### DE STUFF #############################

    def can_run_generations(self):
        return self.it < self.max_iters and self.stop_at_iter

    def can_run_evals(self):
        return self.spent_evals < self.max_evals and self.stop_at_eval

    def population_update(self):
        for i in range(self.pop_size):
            self.target = i
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

    def selection(self, candidate=None):
        self.standard_selection(candidate)

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
        assert self.sade_k is not None, 'sade_k should not be None'
        assert self.target is not None, 'target should not be None'

        return self.sade_cr[self.target][self.sade_k]

    def update_mean(self):
        self.mean = 0
        self.best_score = float('inf')
        for i in range(self.pop_size):
            self.mean += self.pop[i].score / self.pop_size
            if self.best_score is None or self.pop[i].score < self.best_score:
                self.best_score = self.pop[i].score
                self.best_index = i

        return self.mean

    def update_moment_of_inertia(self):
        pop_size = self.pop_size
        n_dim = self.pop[0].total_number_of_angles

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

        total_number_of_angles = len(self.pop[0].angles)

        for i in range(0, self.pop_size):
            for j in range(i + 1, self.pop_size):
                aux_1 = 0

                ind_a = self.pop[i].angles
                ind_b = self.pop[j].angles

                for d in range(0, total_number_of_angles):
                    aux_1 += (ind_a[d] - ind_b[d]) ** 2

                aux_1 = math.sqrt(aux_1) / total_number_of_angles

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

        name = self.rosetta_pack.protein_loader.original + '/' + \
            ("popdata_%08d_" % self.it) + self.name_suffix + ".dat"
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
            s += self.rosetta_pack.get_rmsd_from_native(self.pop[i].pose)

        return s / c

    def log(self, it=None):
        if it is None:
            it = self.it

        rmsd = self.rosetta_pack.get_rmsd_from_native(self.pop[self.best_index].pose)

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

            eta_iters = (self.max_iters - it) * secs_per_iter

            self.last_time = now
            self.last_spent_evals = self.spent_evals
            self.time_pivot += 1

        cr = ''
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

        for k, (a, b) in enumerate(data):
            try:
                string += a % b
                string += ' '
            except Exception:
                print(k, a, b)

        print(string)
        self.stats.write(string + '\n')

        if self.ops_stats is not None and self.sade_success_memory is not None:
            string = ''

            string += '%8d  ' % it

            for k in range(self.sade_n_ops):
                s_s = sum([self.sade_success_memory[i][k] for i in range(self.sade_lp)])
                s_f = sum([self.sade_failure_memory[i][k] for i in range(self.sade_lp)])
                string += '  %8d %8d %8d' % (s_s + s_f, s_s, s_f)
            self.ops_stats.write(string + '\n')
            self.ops_stats.flush()
