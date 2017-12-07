import rosetta_pack
import protein_data
import random
import numpy as np
import math
import sys
import time
import datetime
import string


class DE:
    def __init__(self, pop_size=50, pname='1zdd', c_rate=0.5, f_factor=0.5, max_iters=100, allatom=False):
        self.rosetta_pack = rosetta_pack.RosettaPack(pname)

        self.pname = pname

        self.pop_size = pop_size
        self.pop = [protein_data.ProteinData(self.rosetta_pack, allatom=allatom) for _ in range(pop_size)]
        self.c_rate = c_rate
        self.f_factor = f_factor

        self.trial = protein_data.ProteinData(self.rosetta_pack, allatom=allatom)
        self.max_iters = max_iters
        self.spent_iters = 0
        self.m_nmdf = 0
        self.cname = ''

        self.d = self.pop[0].pose.total_residue()

        # Other stuff
        self.coil_only = False
        self.allatom = False

        self.stage0_init = False
        self.stage0_init = False
        self.stage2_interval = -1
        self.stage2_all_interval = -1
        self.partial_reset = -1
        self.log_interval = 10

        self.reset_d_trigger = 0.0
        self.reset_d_percent = 0.75

        # Island stuff
        self.comm = None  # Comunicator
        self.island_interval = 100

        # LHS parameters
        self.do_lhs = False
        self.n_hashes = 10
        self.hashes = None
        self.hash_values = None
        self.active_hash1 = 0
        self.active_hash2 = 1
        self.tmp1 = None
        self.tmp2 = None
        self.update_interval = 20
        self.change_interval = 100

        # Log stuff
        self.init_time = time.time()
        self.now = datetime.datetime.now()

        # SaDE stuff
        self.sade_run = True
        self.sade_lp = 50
        self.sade_lp_left = self.sade_lp
        self.sade_f = []

        self.sade_ops = [self.best1bin_global, self.rand1bin_global, self.rand2bin_global,
                         self.currToRand_global, self.currToBest_global]
        # self.sade_ops = [self.rand1bin_global, self.rand2bin_global]
        # self.sade_ops = [self.rand1bin_global]
        self.sade_n_ops = len(self.sade_ops)

        self.sade_ops_probs = None #  [1 / self.sade_n_ops for _ in range(self.sade_n_ops)]

        self.sade_success_memory = None #  [[0 for k in range(self.sade_n_ops)] for i in range(self.sade_lp)]
        self.sade_failure_memory = None #  [[0 for k in range(self.sade_n_ops)] for i in range(self.sade_lp)]

        self.sade_cr = None #  [[random.random() for k in range(self.sade_n_ops)] for i in range(self.pop_size)]
        self.sade_cr_m = None #  [random.random() for k in range(self.sade_n_ops)]
        self.sade_cr_memory = None #  [[] for k in range(self.sade_n_ops)]

        self.sade_reinit_interval = 1000

        # Inner info
        self.mean = 0
        self.best_index = 0
        self.best_score = 0

        print('Finished initialization')

    def sade_reinit(self):
        self.sade_cr = [[random.random() for k in range(self.sade_n_ops)] for i in range(self.pop_size)]
        # self.sade_cr_m = [random.random() for k in range(self.sade_n_ops)]  # uniform
        self.sade_cr_m = [np.clip(random.gauss(0.75, 0.1), 0.6, 1.0) for k in range(self.sade_n_ops)]  # gaussian
        self.sade_cr_memory = [[self.sade_cr_m[k]] for k in range(self.sade_n_ops)]

        self.sade_success_memory = [[0 for k in range(self.sade_n_ops)] for i in range(self.sade_lp)]
        self.sade_failure_memory = [[0 for k in range(self.sade_n_ops)] for i in range(self.sade_lp)]

        self.sade_ops_probs = [1 / self.sade_n_ops for _ in range(self.sade_n_ops)]

        self.sade_lp_left = self.sade_lp

    def sade_update_parameters(self):
        self.sade_f = [random.gauss(0.5, 0.3) for _ in range(self.pop_size)]

        if self.sade_lp_left <= 0:
            for k in range(self.sade_n_ops):
                self.sade_cr_memory[k].sort()

            for k in range(self.sade_n_ops):
                if len(self.sade_cr_memory[k]) == 0:
                    print('Empty memory for k = %d %s!' % (k, self.sade_ops[k]))
                else:
                    self.sade_cr_m[k] = np.median(self.sade_cr_memory[k])
                    # if len(self.sade_cr_memory[k]) > 2000:
                    #     m = np.median(self.sade_cr_memory[k])
                    #     self.sade_cr_memory[k] = [m]

            self.sade_cr = [[np.clip(random.gauss(self.sade_cr_m[k], 0.1), 0.0, 1.0) for k in range(self.sade_n_ops)]
                            for i in range(self.pop_size)]
        else:
            self.sade_lp_left -= 1

    def sade_update_ops(self):
        pass

    def sade_get_op(self):
        n = random.random()
        a = 0.0
        i = -1

        while a < n:
            i += 1
            a += self.sade_ops_probs[i]

        return self.sade_ops[i]

    def open_stats(self):
        char_set = string.ascii_uppercase + string.digits
        r_string = ''.join(random.sample(char_set * 6, 6))

        now = self.now

        self.name_suffix = "_%s__%s__%04d_%02d_%02d__%02d_%02d_%02d__%s" % \
                           (self.pname, self.cname, now.year, now.month, now.day, now.hour, now.minute, now.second, r_string)

        self.stats = open(self.rosetta_pack.protein_loader.original + '/' + "stats_" + self.name_suffix + ".dat", 'w')

    def dump_config(self):
        with open(self.rosetta_pack.protein_loader.original + '/' + "parameters_" + self.name_suffix + ".yaml", 'w') as f:
            f.write('pname: %s\n' % (self.pname))
            f.write('pop_size: %d\n' % (self.pop_size))
            f.write('c_rate: %f\n' % (self.c_rate))
            f.write('f_factor: %f\n' % (self.f_factor))
            f.write('max_iters: %d\n' % (self.max_iters))
            f.write('coil_only: %d\n' % (self.coil_only))
            f.write('allatom: %d\n' % (self.allatom))
            f.write('stage0_init: %d\n' % self.stage0_init)
            f.write('stage2_interval: %d\n' % self.stage2_interval)
            f.write('stage2_all_interval: %d\n' % self.stage2_all_interval)
            f.write('partial_reset: %d\n' % self.partial_reset)
            f.write('log_interval: %d\n' % self.log_interval)
            f.write('island_interval: %d\n' % self.island_interval)
            f.write('do_lhs: %d\n' % self.do_lhs)
            f.write('n_hashes: %d\n' % self.n_hashes)
            f.write('update_interval: %d\n' % self.update_interval)
            f.write('change_interval: %d\n' % self.change_interval)
            f.write('reset_d_trigger: %d\n' % self.reset_d_trigger)
            f.write('reset_d_percent: %d\n' % self.reset_d_percent)
            f.write('sade_run: %d\n' % self.sade_run)
            f.write('sade_n_ops: %d\n' % self.sade_n_ops)
            f.write('sade_lp: %d\n' % self.sade_lp)

            f.flush()

    def set_coms(self, pigeon):
        self.comm = pigeon

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

    def create_hashs(self):
        self.hashes = [np.random.randint(100, size=self.pop[0].nsca) for _ in range(self.n_hashes)]

    def apply_hash(self, debug=False):
        # debug = True

        if self.hash_values is None:
            self.hash_values = [[] for _ in range((self.n_hashes + 1) ** 2 + (self.n_hashes + 1))]
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
            tmp1[i] = abs(h1)
            tmp2[i] = abs(h2)

            if minh1 is None or h1 < minh1:
                minh1 = h1

            if maxh1 is None or h1 > maxh1:
                maxh1 = h1

            if minh2 is None or h2 < minh2:
                minh2 = h2

            if maxh2 is None or h2 > maxh2:
                maxh2 = h2

        r1 = (maxh1 - minh1) / self.n_hashes
        b1 = random.random() * r1

        r2 = (maxh2 - minh2) / self.n_hashes
        b2 = random.random() * r2

        if debug:
            print("r: %8.3f  b: %8.3f  minh: %8.3f  maxh: %8.3f" % (r1, b1, minh1, maxh1))
            print("r: %8.3f  b: %8.3f  minh: %8.3f  maxh: %8.3f" % (r2, b2, minh2, maxh2))

        for i in range(self.pop_size):
            a, b = math.floor((tmp1[i] + b1) / r1), math.floor((tmp2[i] + b2) / r2)
            v = (self.n_hashes) * a + b

            if v >= len(self.hash_values):
                print(v, len(self.hash_values))
                v -= len(self.hash_values) - 2

            self.hash_values[v].append(i)

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

        while old2 == self.active_hash2:
            self.active_hash2 = random.randint(0, self.n_hashes - 1)

        if self.active_hash1 == self.active_hash2:
            self.change_hash()

    def get_best(self):
        return self.pop[self.best_index].angles

    def run(self):
        self.open_stats()
        self.dump_config()

        self.create_hashs()

        self.mean = 0
        self.best_index = 0
        self.best_score = None

        for i in range(self.pop_size):
            self.mean += self.pop[i].score / self.pop_size
            if self.best_score is None or self.pop[i].score < self.best_score:
                self.best_score = self.pop[i].score
                self.best_index = i

        self.log(it=-1)

        if self.stage0_init:
            for i in range(self.pop_size):
                if random.random() < .1:
                    self.pop[i].eval()
                    self.pop[i].stage1_mc()
                    self.pop[i].update_angle_from_pose()
                    self.pop[i].eval()

        self.mean = 0
        self.best_index = 0
        self.best_score = None

        for i in range(self.pop_size):
            self.mean += self.pop[i].score / self.pop_size
            if self.best_score is None or self.pop[i].score < self.best_score:
                self.best_score = self.pop[i].score
                self.best_index = i

        self.log(it=0)

        if self.do_lhs:
            self.apply_hash()

        self.start_time = time.time()

        it = 0
        while it < self.max_iters:
            if it % self.sade_reinit_interval == 0 or it == 0:
                self.sade_reinit()

            self.sade_update_parameters()
            self.sade_update_ops()

            it += 1
            self.it = it
            self.best_score = None

            self.mean = 0
            self.best_index = 0

            if self.do_lhs and it % self.change_interval == 0:
                self.change_hash()

            if self.do_lhs and it % self.update_interval == 0:
                self.apply_hash()

            for i in range(self.pop_size):
                if self.do_lhs:
                    self.rand1bin_lhs(i)
                else:
                    if self.sade_run:
                        self.sade_get_op()(i)
                    else:
                        self.rand1bin_global(i)
                # self.mean += self.pop[i].score / self.pop_size

                if self.best_score is None or self.pop[i].score < self.best_score:
                    self.best_score = self.pop[i].score
                    self.best_index = i

            if self.do_lhs and False:
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

            if self.diversity < self.reset_d_trigger:
                print('reset')
                for i in range(self.pop_size):
                    if random.random() < self.reset_percent and i != self.best_index:
                        self.pop[i].reset()
                        self.pop[i].stage1_mc()
                        self.pop[i].update_angle_from_pose()
                        self.pop[i].eval()

            if (self.partial_reset > 0 and it % self.partial_reset == 0 and it > 0):
                print('Partial reset')
                for i in range(self.pop_size):
                    if random.random() < .15 and i != self.best_index:
                        self.pop[i].reset()
                        self.pop[i].stage2_mc()
                        self.pop[i].update_angle_from_pose()
                        self.pop[i].eval()

            if self.stage2_interval > 0 and it % self.stage2_interval == 0 and it > 0:
                print('LS')
                self.pop[self.best_index].stage2_mc()
                self.pop[self.best_index].update_angle_from_pose()
                self.pop[self.best_index].eval()

            if self.stage2_all_interval > 0 and it % self.stage2_all_interval == 0 and it > 0:
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

            for i in range(self.pop_size):
                self.mean += self.pop[i].score / self.pop_size
                if self.best_score is None or self.pop[i].score < self.best_score:
                    self.best_score = self.pop[i].score
                    self.best_index = i

            if self.comm.size > 1 and self.island_interval > 0 and it % self.island_interval == 0 and it > 0:
                print("% is sending obj with score %f" % (self.comm.rank, self.best_score))
                new_guy = self.comm.migration(self.get_best())
                if new_guy is not None:
                    # self.pop[self.best_index].new_angles(new_guy)
                    # self.pop[self.best_index].eval()
                    self.pop[0].new_angles(new_guy)
                    self.pop[0].eval()
                    # print('Ha! migration')

            if it % 100 == 0:
                self.dump_pbd_best(it)

            if self.log_interval > 0 and it % self.log_interval == 0:
                # self.pop[0].print_angles()
                self.log()
                self.stats.flush()
                if self.do_lhs:
                    self.print_hash()

                sys.stdout.flush()

            self.rosetta_pack.pymover.apply(self.pop[self.best_index].pose)

        # end_time = time.time()
        # print("Processing took %f seconds" % (end_time - start_time))

    def print_hash(self):
        for n, i in enumerate(self.hash_values):
            if len(i) > 0:
                print(n, i)

    def rand1bin_lhs(self, huehue):
        hi = 0

        for n, hs in enumerate(self.hash_values):
            if huehue in hs:
                hi = n

        if len(self.hash_values[hi]) < 3:
            return

        ps = random.sample(self.hash_values[hi], k=3)

        p1 = ps[0]
        p2 = ps[1]
        p3 = ps[2]

        cutPoint = random.randint(0, self.rosetta_pack.pose.total_residue())

        t_angle = []

        ind1 = self.pop[p1]
        ind2 = self.pop[p2]
        ind3 = self.pop[p3]

        index = 0
        c = 0
        d = 0
        for k, v in enumerate(self.rosetta_pack.target):
            na = 3 + self.rosetta_pack.bounds.getNumSideChainAngles(v)
            for j in range(na):
                d = index + j
                r = random.random()
                if r < self.c_rate or d == cutPoint:
                    if self.coil_only and self.rosetta_pack.ss_pred[c // 3] != 'C':
                        t_angle.append(ind1.angles[d] + (self.f_factor * (ind2.angles[d] - ind3.angles[d])))
                    elif not self.coil_only:
                        t_angle.append(ind1.angles[d] + (self.f_factor * (ind2.angles[d] - ind3.angles[d])))
                else:
                    t_angle.append(self.pop[huehue].angles[d])

            c += 1
            index += na

        self.trial.new_angles(t_angle)
        self.trial.fix_bounds()
        self.trial.eval()

        if self.trial.score < self.pop[huehue].score:
            t = self.pop[huehue]
            self.pop[huehue] = self.trial
            self.trial = t
            if self.trial is self.pop[huehue]:
                import sys
                sys.exit()

# ########### Global operators

    def best1bin_global(self, huehue):
        sade_k = self.sade_ops.index(self.best1bin_global)

        p1 = self.best_index
        p2 = random.randint(0, self.pop_size - 1)
        p3 = random.randint(0, self.pop_size - 1)

        while p1 == p2 or p2 == p3 or p1 == p3 or p2 == huehue or p3 == huehue:
            p2 = random.randint(0, self.pop_size - 1)
            p3 = random.randint(0, self.pop_size - 1)

        cutPoint = random.randint(0, self.rosetta_pack.pose.total_residue())

        t_angle = []

        ind1 = self.pop[p1]
        ind2 = self.pop[p2]
        ind3 = self.pop[p3]

        index = 0
        c = 0
        d = 0

        f = self.f_factor
        cr = self.c_rate

        if self.sade_run:
            f = self.sade_f[huehue]
            cr = self.sade_cr[huehue][sade_k]

        for k, v in enumerate(self.rosetta_pack.target):
            na = 3 + self.rosetta_pack.bounds.getNumSideChainAngles(v)
            for j in range(na):
                d = index + j
                r = random.random()
                if r < cr or d == cutPoint:
                    if self.coil_only and self.rosetta_pack.ss_pred[c // 3] != 'C':
                        t_angle.append(ind1.angles[d] + (f * (ind2.angles[d] - ind3.angles[d])))
                    elif not self.coil_only:
                        t_angle.append(ind1.angles[d] + (f * (ind2.angles[d] - ind3.angles[d])))
                    else:
                        t_angle.append(self.pop[huehue].angles[d])
                else:
                    t_angle.append(self.pop[huehue].angles[d])

            c += 1
            index += na

        self.trial.new_angles(t_angle)
        self.trial.fix_bounds()
        self.trial.eval()

        if self.trial.score < self.pop[huehue].score:
            if self.sade_run:
                self.sade_cr_memory[sade_k].append(cr)
                ind = self.it % self.sade_lp
                self.sade_success_memory[ind][sade_k] += 1
            t = self.pop[huehue]
            self.pop[huehue] = self.trial
            self.trial = t
            if self.trial is self.pop[huehue]:
                import sys
                sys.exit()
        else:
            if self.sade_run:
                ind = self.it % self.sade_lp
                self.sade_failure_memory[ind][sade_k] += 1

    def rand1bin_global(self, huehue):
        sade_k = self.sade_ops.index(self.rand1bin_global)

        p1 = random.randint(0, self.pop_size - 1)
        p2 = random.randint(0, self.pop_size - 1)
        p3 = random.randint(0, self.pop_size - 1)

        while p1 == p2 or p2 == p3 or p1 == p3 or p1 == huehue or p2 == huehue or p3 == huehue:
            p1 = random.randint(0, self.pop_size - 1)
            p2 = random.randint(0, self.pop_size - 1)
            p3 = random.randint(0, self.pop_size - 1)

        cutPoint = random.randint(0, self.rosetta_pack.pose.total_residue())

        t_angle = []

        ind1 = self.pop[p1]
        ind2 = self.pop[p2]
        ind3 = self.pop[p3]

        index = 0
        c = 0
        d = 0

        f = self.f_factor
        cr = self.c_rate

        if self.sade_run:
            f = self.sade_f[huehue]
            cr = self.sade_cr[huehue][sade_k]

        for k, v in enumerate(self.rosetta_pack.target):
            na = 3 + self.rosetta_pack.bounds.getNumSideChainAngles(v)
            for j in range(na):
                d = index + j
                r = random.random()
                if r < cr:  # or d == cutPoint:
                    if self.coil_only and self.rosetta_pack.ss_pred[c // 3] != 'C':
                        t_angle.append(ind1.angles[d] + (f * (ind2.angles[d] - ind3.angles[d])))
                    elif not self.coil_only:
                        t_angle.append(ind1.angles[d] + (f * (ind2.angles[d] - ind3.angles[d])))
                    else:
                        t_angle.append(self.pop[huehue].angles[d])
                else:
                    t_angle.append(self.pop[huehue].angles[d])

            c += 1
            index += na

        self.trial.new_angles(t_angle)
        self.trial.fix_bounds()
        self.trial.eval()

        if self.trial.score < self.pop[huehue].score:
            if self.sade_run:
                self.sade_cr_memory[sade_k].append(cr)
                ind = self.it % self.sade_lp
                self.sade_success_memory[ind][sade_k] += 1
            t = self.pop[huehue]
            self.pop[huehue] = self.trial
            self.trial = t
            if self.trial is self.pop[huehue]:
                import sys
                sys.exit()
        else:
            if self.sade_run:
                ind = self.it % self.sade_lp
                self.sade_failure_memory[ind][sade_k] += 1

    def rand2bin_global(self, huehue):
        sade_k = self.sade_ops.index(self.rand2bin_global)

        p1 = random.randint(0, self.pop_size - 1)
        p2 = random.randint(0, self.pop_size - 1)
        p3 = random.randint(0, self.pop_size - 1)
        p4 = random.randint(0, self.pop_size - 1)
        p5 = random.randint(0, self.pop_size - 1)

        while p1 == p2 or p2 == p3 or p1 == p3 or p1 == huehue or p2 == huehue or p3 == huehue or p4 == p5 or p4 == huehue or p5 == huehue:
            p1 = random.randint(0, self.pop_size - 1)
            p2 = random.randint(0, self.pop_size - 1)
            p3 = random.randint(0, self.pop_size - 1)
            p4 = random.randint(0, self.pop_size - 1)
            p5 = random.randint(0, self.pop_size - 1)

        cutPoint = random.randint(0, self.rosetta_pack.pose.total_residue())

        t_angle = []

        ind1 = self.pop[p1]
        ind2 = self.pop[p2]
        ind3 = self.pop[p3]
        ind4 = self.pop[p4]
        ind5 = self.pop[p5]

        index = 0
        c = 0
        d = 0

        f = self.f_factor
        cr = self.c_rate

        if self.sade_run:
            f = self.sade_f[huehue]
            cr = self.sade_cr[huehue][sade_k]

        for k, v in enumerate(self.rosetta_pack.target):
            na = 3 + self.rosetta_pack.bounds.getNumSideChainAngles(v)
            for j in range(na):
                d = index + j
                r = random.random()
                if r < cr or d == cutPoint:
                    if self.coil_only and self.rosetta_pack.ss_pred[c // 3] != 'C':
                        t_angle.append(ind1.angles[d] + (f * (ind2.angles[d] - ind3.angles[d])) + (f * (ind4.angles[d] - ind5.angles[d])))
                    elif not self.coil_only:
                        t_angle.append(ind1.angles[d] + (f * (ind2.angles[d] - ind3.angles[d])) + (f * (ind4.angles[d] - ind5.angles[d])))
                    else:
                        t_angle.append(self.pop[huehue].angles[d])
                else:
                    t_angle.append(self.pop[huehue].angles[d])

            c += 1
            index += na

        self.trial.new_angles(t_angle)
        self.trial.fix_bounds()
        self.trial.eval()

        if self.trial.score < self.pop[huehue].score:
            if self.sade_run:
                self.sade_cr_memory[sade_k].append(cr)
                ind = self.it % self.sade_lp
                self.sade_success_memory[ind][sade_k] += 1
            t = self.pop[huehue]
            self.pop[huehue] = self.trial
            self.trial = t
            if self.trial is self.pop[huehue]:
                import sys
                sys.exit()
        else:
            if self.sade_run:
                ind = self.it % self.sade_lp
                self.sade_failure_memory[ind][sade_k] += 1

    def currToRand_global(self, huehue):
        sade_k = self.sade_ops.index(self.currToRand_global)
        sade_k = 3

        p1 = huehue
        p2 = random.randint(0, self.pop_size - 1)
        p3 = random.randint(0, self.pop_size - 1)

        while p1 == p2 or p2 == p3 or p1 == p3 or p2 == huehue or p3 == huehue:
            p1 = random.randint(0, self.pop_size - 1)
            p2 = random.randint(0, self.pop_size - 1)
            p3 = random.randint(0, self.pop_size - 1)

        cutPoint = random.randint(0, self.rosetta_pack.pose.total_residue())

        t_angle = []

        ind1 = self.pop[p1]
        ind2 = self.pop[p2]
        ind3 = self.pop[p3]

        index = 0
        c = 0
        d = 0

        f = self.f_factor
        cr = self.c_rate

        if self.sade_run:
            f = self.sade_f[huehue]
            cr = self.sade_cr[huehue][sade_k]

        for k, v in enumerate(self.rosetta_pack.target):
            na = 3 + self.rosetta_pack.bounds.getNumSideChainAngles(v)
            for j in range(na):
                d = index + j
                r = random.random()
                if r < cr or d == cutPoint:
                    if self.coil_only and self.rosetta_pack.ss_pred[c // 3] != 'C':
                        t_angle.append(ind1.angles[d] + (f * (ind2.angles[d] - ind3.angles[d])))
                    elif not self.coil_only:
                        t_angle.append(ind1.angles[d] + (f * (ind2.angles[d] - ind3.angles[d])))
                    else:
                        t_angle.append(self.pop[huehue].angles[d])
                else:
                    t_angle.append(self.pop[huehue].angles[d])

            c += 1
            index += na

        self.trial.new_angles(t_angle)
        self.trial.fix_bounds()
        self.trial.eval()

        if self.trial.score < self.pop[huehue].score:
            if self.sade_run:
                self.sade_cr_memory[sade_k].append(cr)
                ind = self.it % self.sade_lp
                self.sade_success_memory[ind][sade_k] += 1
            t = self.pop[huehue]
            self.pop[huehue] = self.trial
            self.trial = t
            if self.trial is self.pop[huehue]:
                import sys
                sys.exit()
        else:
            if self.sade_run:
                ind = self.it % self.sade_lp
                self.sade_failure_memory[ind][sade_k] += 1

    def currToBest_global(self, huehue):
        sade_k = self.sade_ops.index(self.currToBest_global)

        p1 = huehue
        p2 = random.randint(0, self.pop_size - 1)
        p3 = random.randint(0, self.pop_size - 1)
        p4 = self.best_index
        p5 = random.randint(0, self.pop_size - 1)

        while p1 == p2 or p2 == p3 or p1 == p3 or p2 == huehue or p3 == huehue or p4 == p5:
            p1 = random.randint(0, self.pop_size - 1)
            p2 = random.randint(0, self.pop_size - 1)
            p3 = random.randint(0, self.pop_size - 1)
            p5 = random.randint(0, self.pop_size - 1)

        cutPoint = random.randint(0, self.rosetta_pack.pose.total_residue())

        t_angle = []

        ind1 = self.pop[p1]
        ind2 = self.pop[p2]
        ind3 = self.pop[p3]
        ind4 = self.pop[p4]
        ind5 = self.pop[p5]

        index = 0
        c = 0
        d = 0

        f = self.f_factor
        cr = self.c_rate

        if self.sade_run:
            f = self.sade_f[huehue]
            cr = self.sade_cr[huehue][sade_k]

        for k, v in enumerate(self.rosetta_pack.target):
            na = 3 + self.rosetta_pack.bounds.getNumSideChainAngles(v)
            for j in range(na):
                d = index + j
                r = random.random()
                if r < cr or d == cutPoint:
                    if self.coil_only and self.rosetta_pack.ss_pred[c // 3] != 'C':
                        t_angle.append(ind1.angles[d] + (f * (ind2.angles[d] - ind3.angles[d])) + (f * (ind4.angles[d] - ind5.angles[d])))
                    elif not self.coil_only:
                        t_angle.append(ind1.angles[d] + (f * (ind2.angles[d] - ind3.angles[d])) + (f * (ind4.angles[d] - ind5.angles[d])))
                    else:
                        t_angle.append(self.pop[huehue].angles[d])
                else:
                    t_angle.append(self.pop[huehue].angles[d])

            c += 1
            index += na

        self.trial.new_angles(t_angle)
        self.trial.fix_bounds()
        self.trial.eval()

        if self.trial.score < self.pop[huehue].score:
            if self.sade_run:
                self.sade_cr_memory[sade_k].append(cr)
                ind = self.it % self.sade_lp
                self.sade_success_memory[ind][sade_k] += 1
            t = self.pop[huehue]
            self.pop[huehue] = self.trial
            self.trial = t
            if self.trial is self.pop[huehue]:
                import sys
                sys.exit()
        else:
            if self.sade_run:
                ind = self.it % self.sade_lp
                self.sade_failure_memory[ind][sade_k] += 1

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

    def dump_pbd_pop(self):
        pass

    def dump_pbd_best(self, it):
        name = self.rosetta_pack.protein_loader.original + '/' + ("best_%05d_" % it) + self.name_suffix + ".pdb"
        self.pop[self.best_index].pose.dump_pdb(name)

    def avg_distance(self):
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
        s = 0
        c = 0

        for i in range(self.pop_size):
            for j in range(i + 1, self.pop_size):
                c += 1
                s += self.rosetta_pack.get_rmsd_from_pose(self.pop[i].pose, self.pop[j].pose)

        return s / c

    def avg_rmsd_from_native(self):
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
            eta = 0
        else:
            secs_per_iter = (time.time() - self.start_time) / it
            eta = (self.max_iters - it) * secs_per_iter

        cr = '%3.2f' % self.c_rate

        if self.sade_run:
            cr = ''

            if self.sade_cr_m is not None:
                for i in range(len(self.sade_cr_m)):
                    cr += '%3.2f ' % self.sade_cr_m[i]

        string = "%2d %8d %12.4f %12.4f %8.4f %8.4f %8.4f %8.4f %7.3f %8.3f %s"
        data = (self.comm.rank, it, self.best_score, self.mean, self.update_diversity(), self.avg_rmsd(),
                rmsd, self.avg_rmsd_from_native(), secs_per_iter, eta, cr)

        # print(self.sade_success_memory)
        # print(self.sade_failure_memory)

        self.stats.write((string + '\n') % data)
        # if it % 100 == 0:
        print(string % data)
