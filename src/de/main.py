from mpi4py import MPI
import rosetta_pack
import protein_data
import random
import sys
import slave
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

        self.comm = MPI.COMM_WORLD
        self.size = comm.size
        self.rank = comm.rank
        self.status = MPI.Status()

    def run(self):
        start_time = time.time()

        self.comm.bcast((self.max_iters, self.c_rate, self.f_factor, self.pop_size, self.rosetta_pack.pose.total_residue(),
                        self.pname))
        for it in range(self.max_iters):
            # print("MASTER: loop")
            newpop = []
            mean = 0
            while len(newpop) < self.pop_size:
                tag, ind = self.comm.recv(source=MPI.ANY_SOURCE, status=self.status)

                if tag == 'BOOT':
                    # print('MASTER got BOOT from %2d' % self.status.Get_source())
                    pass
                else:
                    # print('MASTER got data from %2d with score %8.2f' % (self.status.Get_source(), ind.score))
                    newpop.append(ind)

                    if self.best_score is None or ind.score < self.best_score:
                        self.best_score = ind.score
                    mean += ind.score / self.pop_size

                data = None

                p1 = random.randint(0, self.pop_size - 1)
                p2 = random.randint(0, self.pop_size - 1)
                p3 = random.randint(0, self.pop_size - 1)

                while p1 == p2 or p2 == p3 or p1 == p3:
                    p1 = random.randint(0, self.pop_size - 1)
                    p2 = random.randint(0, self.pop_size - 1)
                    p3 = random.randint(0, self.pop_size - 1)

                cutPoint = random.randint(0, self.rosetta_pack.pose.total_residue())

                ind1 = self.pop[p1]
                ind2 = self.pop[p2]
                ind3 = self.pop[p3]

                data = ('WORK', ind1, ind2, ind3, cutPoint)

                self.comm.send(data, dest=self.status.Get_source())

            self.pop = newpop
            print("%8d %8.3f %8.3f" % (it, self.best_score, mean))
            sys.stdout.flush()

        for _ in range(self.size - 1):
            tag, ind = self.comm.recv(source=MPI.ANY_SOURCE, status=self.status)
            self.comm.send(('DIE', None, None, None, None), dest=self.status.Get_source())

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
    comm = MPI.COMM_WORLD
    size = comm.size
    rank = comm.rank

    if rank == 0:
        de = DE(c_rate=0.1, f_factor=0.9, max_iters=10)
        de.run()
    else:
        slave.run()
