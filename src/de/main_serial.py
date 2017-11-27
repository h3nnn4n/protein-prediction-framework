import rosetta_pack
import protein_data
import random
import sys
import time


class DE:
    def __init__(self, pop_size=50, pname='1zdd', c_rate=0.9, f_factor=0.8, max_iters=100):
        self.rosetta_pack = rosetta_pack.RosettaPack(pname)

        self.pname = pname

        self.pop_size = pop_size
        self.pop = [protein_data.ProteinData(self.rosetta_pack) for _ in range(pop_size)]
        self.c_rate = c_rate
        self.f_factor = f_factor

        self.max_iters = max_iters

        self.best_score = None

    def run(self):
        start_time = time.time()

        for it in range(self.max_iters):
            newpop = []
            self.best_score = None

            mean = 0

            for i in range(self.pop_size):
                trial = self.rand1bin()
                trial.eval()

                if trial.score < self.pop[i].score:
                    # print("%8.3f %8.3f -- NEW" % (trial.score, self.pop[i].score))
                    newpop.append(trial)
                    self.rosetta_pack.pymover.apply(trial.pose)
                else:
                    # print("%8.3f %8.3f -- SAME" % (trial.score, self.pop[i].score))
                    newpop.append(self.pop[i])

                if self.best_score is None or newpop[i].score < self.best_score:
                    self.best_score = newpop[i].score
                    # self.rosetta_pack.pymover.apply(newpop[i].pose)
                    # print("         %8.3f" % (self.best_score))

                mean += newpop[i].score / self.pop_size

            print("%8d %8.3f %8.3f" % (it, self.best_score, mean))
            sys.stdout.flush()

        end_time = time.time()
        print("Processing took %f seconds" % (end_time - start_time))

    def rand1bin(self):
        p1 = random.randint(0, self.pop_size - 1)
        p2 = random.randint(0, self.pop_size - 1)
        p3 = random.randint(0, self.pop_size - 1)

        while p1 == p2 or p2 == p3 or p1 == p3:
            p1 = random.randint(0, self.pop_size - 1)
            p2 = random.randint(0, self.pop_size - 1)
            p3 = random.randint(0, self.pop_size - 1)

        cutPoint = random.randint(0, self.rosetta_pack.pose.total_residue())

        t_phi = []
        t_psi = []
        t_omega = []

        ind1 = self.pop[p1]
        ind2 = self.pop[p2]
        ind3 = self.pop[p3]

        for d in range(self.rosetta_pack.pose.total_residue() - 1):
            if random.random() < self.c_rate or d == cutPoint:
                t_phi.append(ind1.phi[d] + (self.f_factor * (ind2.phi[d] - ind3.phi[d])))
                t_psi.append(ind1.psi[d] + (self.f_factor * (ind2.psi[d] - ind3.psi[d])))
                t_omega.append(ind1.omega[d] + (self.f_factor * (ind2.omega[d] - ind3.omega[d])))
            else:
                t_phi.append(ind1.phi[d])
                t_psi.append(ind1.psi[d])
                t_omega.append(ind1.omega[d])

        trial = protein_data.ProteinData(self.rosetta_pack)
        trial.new_angles(t_phi, t_psi, t_omega)
        trial.fix_bounds()
        return trial


if __name__ == '__main__':
    de = DE(c_rate=0.1, f_factor=0.9, max_iters=10)
    de.run()
