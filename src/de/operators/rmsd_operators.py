import random
from operators.global_exp import exp_base_operator


def rand1exp_global_max_rmsd(self, target):
    individuals = self.get_individual_within_rmsd_range(n=3, min_range=2, max_range=5)

    if len(individuals) != 3:
        return 0

    exp_base_operator(self, individuals, target, self.rand1exp_global_max_rmsd)


def best1exp_global_max_rmsd(self, target):
    best = self.get_best_individual()
    individuals = self.get_individual_within_rmsd_range(n=3, min_range=2, max_range=5, first=best)

    if len(individuals) != 3:
        return 0

    exp_base_operator(self, individuals, target, self.best1exp_global_max_rmsd)
