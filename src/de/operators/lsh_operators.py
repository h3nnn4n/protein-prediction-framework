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
    # self.de.trial.fix_bounds()
    self.de.trial.eval()

    self.de.selection(self.de.pop[target])


def rand1exp_lsh(self, target):
    if not self.de.sade_run:
        sade_k = 0
    else:
        sade_k = self.de.sade_ops.index(self.rand1exp_lsh)
    self.de.sade_k = sade_k
    hi = 0

    if self.de.locality_sensitive_hashing.hash_values is None:
        raise Exception('Attempted to use rand1exp_lsh without initializing lsh!')

    for n, hs in enumerate(self.de.locality_sensitive_hashing.hash_values):
        if target in hs:
            hi = n

    if len(self.de.locality_sensitive_hashing.hash_values[hi]) < 3:
        return

    ps = random.sample(self.de.locality_sensitive_hashing.hash_values[hi], k=3)

    p1 = ps[0]
    p2 = ps[1]
    p3 = ps[2]

    cutPoint = random.randint(0, self.de.pop[0].total_number_of_angles)

    t_angle = []

    ind1 = self.de.pop[p1]
    ind2 = self.de.pop[p2]
    ind3 = self.de.pop[p3]

    f, cr = self.de.get_f_cr()

    L = 0
    r = 0.0
    pivot = cutPoint

    for i in range(0, ind1.total_number_of_angles):
        t_angle.append(self.de.pop[target].angles[i])

    while L < ind1.total_number_of_angles and r < cr:
        t_angle[pivot % ind1.total_number_of_angles] = (
            ind1.angles[pivot % ind1.total_number_of_angles] +
            (f * (ind2.angles[pivot % ind1.total_number_of_angles] - ind3.angles[pivot % ind1.total_number_of_angles]))
        )

        r = random.random()
        L += 1
        pivot += 1

    self.de.trial.new_angles(t_angle)
    # self.de.trial.fix_bounds()
    self.de.trial.eval()

    self.de.selection(self.de.pop[target])
