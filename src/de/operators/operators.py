import random
import types
from operators.rand1bin_global import rand1bin_global
from operators.rand2bin_global import rand2bin_global
from operators.best1bin_global import best1bin_global
from operators.best2bin_global import best2bin_global
from operators.currToBest_global import currToBest_global
from operators.currToRand_global import currToRand_global
from operators.global_exp import currToRand_exp_global, currToBest_exp_global, \
    rand1exp_global, rand2exp_global, best1exp_global, best2exp_global
from operators.rosetta_operators import monte_carlo_3, monte_carlo_3s, \
    monte_carlo_9, monte_carlo_9s, monte_carlo_3_coil_only
from operators.lsh_operators import rand1exp_lsh
from operators.rmsd_operators import rand1exp_global_max_rmsd, \
    best1exp_global_max_rmsd
from operators.hooke_jeeves_operators import hooke_jeeves_ls


class Operators:
    def __init__(self, de=None):
        self.de = de
        self.available_operators = {
            'rand1bin_global': self.get_rand1bin_global,
            'rand2bin_global': self.get_rand2bin_global,
            'best1bin_global': self.get_best1bin_global,
            'best2bin_global': self.get_best2bin_global,
            'currToBest_global': self.get_currToBest_global,
            'currToRand_global': self.get_currToRand_global,
            'currToRand_exp_global': self.get_currToRand_exp_global,
            'currToBest_exp_global': self.get_currToBest_exp_global,
            'rand1exp_global': self.get_rand1exp_global,
            'rand2exp_global': self.get_rand2exp_global,
            'best1exp_global': self.get_best1exp_global,
            'best2exp_global': self.get_best2exp_global,
            'monte_carlo_3': self.get_monte_carlo_3,
            'monte_carlo_3s': self.get_monte_carlo_3s,
            'monte_carlo_9': self.get_monte_carlo_9,
            'monte_carlo_9s': self.get_monte_carlo_9s,
            'monte_carlo_3_coil_only': self.get_monte_carlo_3_coil_only,
            'rand1exp_lsh': self.get_rand1exp_lsh,
            'rand1exp_global_max_rmsd': self.get_rand1exp_global_max_rmsd,
            'best1exp_global_max_rmsd': self.get_best1exp_global_max_rmsd,
            'hooke_jeeves_ls': self.get_hooke_jeeves_ls,
        }

    def get_random_individual(self):
        p = random.randint(0, self.de.pop_size - 1)
        return self.de.pop[p]

    def get_random_individuals(self, n=1):
        return random.sample(self.de.pop, n)

    def get_best_individual(self):
        return self.de.pop[self.de.best_index]

    def get_current_individual(self, current):
        return self.de.pop[current]

    def get_individual_within_rmsd_range(self, n=3, min_range=1, max_range=50, n_retries=10, first=None):
        for _ in range(n_retries):
            if first is None:
                first = self.get_random_individual()

            individuals_list = [first]

            for _ in range(1, n):
                for _, candidate in enumerate(self.de.pop):
                    rmsd = self.de.rosetta_pack.get_rmsd_from_pose(first.pose, candidate.pose)
                    if min_range <= rmsd and rmsd <= max_range and candidate not in individuals_list:
                        individuals_list.append(candidate)

                    if len(individuals_list) == n:
                        return individuals_list

        return []

    def get_operator(self, name):
        return self.available_operators[name]()

    def get_rand1bin_global(self):
        self.rand1bin_global = types.MethodType(rand1bin_global, self)
        return self.rand1bin_global

    def get_rand2bin_global(self):
        self.rand2bin_global = types.MethodType(rand2bin_global, self)
        return self.rand2bin_global

    def get_best1bin_global(self):
        self.best1bin_global = types.MethodType(best1bin_global, self)
        return self.best1bin_global

    def get_best2bin_global(self):
        self.best2bin_global = types.MethodType(best2bin_global, self)
        return self.best2bin_global

    def get_currToBest_global(self):
        self.currToBest_global = types.MethodType(currToBest_global, self)
        return self.currToBest_global

    def get_currToRand_global(self):
        self.currToRand_global = types.MethodType(currToRand_global, self)
        return self.currToRand_global

    def get_currToRand_exp_global(self):
        self.currToRand_exp_global = types.MethodType(currToRand_exp_global, self)
        return self.currToRand_exp_global

    def get_currToBest_exp_global(self):
        self.currToBest_exp_global = types.MethodType(currToBest_exp_global, self)
        return self.currToBest_exp_global

    def get_rand1exp_global(self):
        self.rand1exp_global = types.MethodType(rand1exp_global, self)
        return self.rand1exp_global

    def get_rand2exp_global(self):
        self.rand2exp_global = types.MethodType(rand2exp_global, self)
        return self.rand2exp_global

    def get_best1exp_global(self):
        self.best1exp_global = types.MethodType(best1exp_global, self)
        return self.best1exp_global

    def get_best2exp_global(self):
        self.best2exp_global = types.MethodType(best2exp_global, self)
        return self.best2exp_global

    def get_monte_carlo_3(self):
        self.monte_carlo_3 = types.MethodType(monte_carlo_3, self)
        return self.monte_carlo_3

    def get_monte_carlo_3s(self):
        self.monte_carlo_3s = types.MethodType(monte_carlo_3s, self)
        return self.monte_carlo_3s

    def get_monte_carlo_9(self):
        self.monte_carlo_9 = types.MethodType(monte_carlo_9, self)
        return self.monte_carlo_9

    def get_monte_carlo_9s(self):
        self.monte_carlo_9s = types.MethodType(monte_carlo_9s, self)
        return self.monte_carlo_9s

    def get_monte_carlo_3_coil_only(self):
        self.monte_carlo_3_coil_only = types.MethodType(monte_carlo_3_coil_only, self)
        return self.monte_carlo_3_coil_only

    def get_rand1exp_lsh(self):
        self.rand1exp_lsh = types.MethodType(rand1exp_lsh, self)
        return self.rand1exp_lsh

    def get_rand1exp_global_max_rmsd(self):
        self.rand1exp_global_max_rmsd = types.MethodType(rand1exp_global_max_rmsd, self)
        return self.rand1exp_global_max_rmsd

    def get_best1exp_global_max_rmsd(self):
        self.best1exp_global_max_rmsd = types.MethodType(best1exp_global_max_rmsd, self)
        return self.best1exp_global_max_rmsd

    def get_hooke_jeeves_ls(self):
        self.hooke_jeeves_ls = types.MethodType(hooke_jeeves_ls, self)
        return self.hooke_jeeves_ls
