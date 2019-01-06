import random


def currToRand_exp_global(self, target):
    if not self.de.sade_run:
        sade_k = 0
    else:
        sade_k = self.de.sade_ops.index(self.currToRand_exp_global)
    self.de.sade_k = sade_k

    cutPoint = random.randint(0, self.de.pop[0].nsca)
    t_angle = []
    ind1 = self.get_current_individual(target)
    ind2, ind3 = self.get_random_individuals(n=2)
    f, cr = self.de.get_f_cr()
    L, r, pivot = 0, 0.0, cutPoint

    for i in range(0, ind1.nsca):
        t_angle.append(self.de.pop[target].angles[i])

    while L < ind1.nsca and r < cr:
            t_angle[pivot % ind1.nsca] = (
                ind1.angles[pivot % ind1.nsca] +
                (f * (ind2.angles[pivot % ind1.nsca] - ind3.angles[pivot % ind1.nsca]))
            )

            r = random.random()
            L += 1
            pivot += 1

    self.de.trial.new_angles(t_angle)
    self.de.trial.fix_bounds()
    self.de.trial.eval()

    self.de.selection(self.de.pop[target])


def currToBest_exp_global(self, target):
    if not self.de.sade_run:
        sade_k = 0
    else:
        sade_k = self.de.sade_ops.index(self.currToBest_exp_global)
    self.de.sade_k = sade_k

    cutPoint = random.randint(0, self.de.pop[0].nsca)
    t_angle = []
    ind1 = self.get_current_individual(target)
    ind2 = self.get_best_individual()
    ind3 = self.get_random_individual()
    f, cr = self.de.get_f_cr()
    L, r, pivot = 0, 0.0, cutPoint

    for i in range(0, ind1.nsca):
        t_angle.append(self.de.pop[target].angles[i])

    while L < ind1.nsca and r < cr:
            t_angle[pivot % ind1.nsca] = (
                ind1.angles[pivot % ind1.nsca] +
                (f * (ind2.angles[pivot % ind1.nsca] - ind3.angles[pivot % ind1.nsca]))
            )

            r = random.random()
            L += 1
            pivot += 1

    self.de.trial.new_angles(t_angle)
    self.de.trial.fix_bounds()
    self.de.trial.eval()

    self.de.selection(self.de.pop[target])
