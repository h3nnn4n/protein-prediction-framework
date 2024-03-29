import random


def exp_base_operator(self, individuals, target, operator):
    if not self.de.sade_run:
        sade_k = 0
    else:
        sade_k = self.de.sade_ops.index(operator)
    self.de.sade_k = sade_k

    cutPoint = random.randint(0, self.de.pop[0].total_number_of_angles)
    t_angle = []
    f, cr = self.de.get_f_cr()
    L, r, pivot = 0, 0.0, cutPoint

    for i in range(0, individuals[0].total_number_of_angles):
        t_angle.append(self.de.pop[target].angles[i])

    if len(individuals) == 3:
        ind1, ind2, ind3 = individuals
        while L < ind1.total_number_of_angles and r < cr:
            index = pivot % ind1.total_number_of_angles
            t_angle[index] = (
                ind1.angles[index] +
                (f * (ind2.angles[index] - ind3.angles[index]))
            )

            r = random.random()
            L += 1
            pivot += 1

    elif len(individuals) == 5:
        ind1, ind2, ind3, ind4, ind5 = individuals
        f1, f2 = self.de.get_f(), self.de.get_f()
        while L < ind1.total_number_of_angles and r < cr:
            index = pivot % ind1.total_number_of_angles
            t_angle[index] = (
                ind1.angles[index] +
                (f1 * (ind2.angles[index] - ind3.angles[index])) +
                (f2 * (ind4.angles[index] - ind5.angles[index]))
            )

            r = random.random()
            L += 1
            pivot += 1

    else:
        raise ValueError('Invalid number of individuals, please use 3 or 5')

    self.de.trial.new_angles(t_angle)
    self.de.trial.fix_bounds()
    self.de.trial.eval()

    self.de.selection(self.de.pop[target])


def rand1exp_global(self, target):
    individuals = self.get_random_individuals(n=3)
    exp_base_operator(self, individuals, target, self.rand1exp_global)


def rand2exp_global(self, target):
    individuals = self.get_random_individuals(n=5)
    exp_base_operator(self, individuals, target, self.rand2exp_global)


def best1exp_global(self, target):
    individuals = [self.get_best_individual()]
    individuals.extend(self.get_random_individuals(n=2))
    exp_base_operator(self, individuals, target, self.best1exp_global)


def best2exp_global(self, target):
    individuals = [self.get_best_individual()]
    individuals.extend(self.get_random_individuals(n=4))
    exp_base_operator(self, individuals, target, self.best2exp_global)


def currToRand_exp_global(self, target):
    individuals = [self.get_current_individual(target)]
    individuals.extend(self.get_random_individuals(n=2))
    exp_base_operator(self, individuals, target, self.currToRand_exp_global)


def currToBest_exp_global(self, target):
    individuals = [
        self.get_current_individual(target),
        self.get_best_individual(),
        self.get_random_individual()
    ]
    # individuals.extend(self.get_random_individuals(n=2))
    exp_base_operator(self, individuals, target, self.currToBest_exp_global)
