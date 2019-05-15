import random
from local_search.hooke_jeeves import hooke


def hooke_jeeves_ls(self, target):
    if not self.de.sade_run:
        sade_k = 0
    else:
        sade_k = self.de.sade_ops.index(self.hooke_jeeves_ls)
    self.de.sade_k = sade_k

    individual = self.get_current_individual(target)

    f, _ = self.de.get_f_cr()

    evals, _ = hooke(individual, rho=f)

    self.de.selection(self.de.pop[target])

    return evals
