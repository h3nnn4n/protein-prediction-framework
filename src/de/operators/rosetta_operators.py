def base_monte_carlo(self, target, mode, operator):
    if not self.de.sade_run:
        sade_k = 0
    else:
        sade_k = self.de.sade_ops.index(operator)
    self.de.sade_k = sade_k

    _, cr = self.de.get_f_cr()
    self.de.trial.copy(self.de.pop[target])
    evals = self.de.trial.stage2_mc(n=25, temp=cr * 3.0, mode=mode)

    self.de.selection(self.de.pop[target])

    return evals


def monte_carlo_3(self, target):
    return base_monte_carlo(self, target, '3', self.monte_carlo_3)


def monte_carlo_3s(self, target):
    return base_monte_carlo(self, target, '3s', self.monte_carlo_3s)


def monte_carlo_9(self, target):
    return base_monte_carlo(self, target, '9', self.monte_carlo_9)


def monte_carlo_9s(self, target):
    return base_monte_carlo(self, target, '9s', self.monte_carlo_9s)


def base_coil_search(self, target, mode, operator):
    if not self.de.sade_run:
        sade_k = 0
    else:
        sade_k = self.de.sade_ops.index(operator)
    self.de.sade_k = sade_k

    old_mm = self.de.rosetta_pack.mover_3mer.movemap()

    _, cr = self.de.get_f_cr()
    self.de.trial.copy(self.de.pop[target])
    evals = self.de.trial.stage2_mc(n=10, temp=cr * 3.0, mode=mode)
    self.de.rosetta_pack.mover_3mer.set_movemap(old_mm)

    self.de.selection(self.de.pop[target])

    return evals


def monte_carlo_3_coil_only(self, target):
    return base_coil_search(self, target, '3', self.monte_carlo_3_coil_only)
