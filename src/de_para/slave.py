from mpi4py import MPI
import random
import protein_data
import rosetta_pack


comm = MPI.COMM_WORLD
size = comm.size
rank = comm.rank


def run():
    # print("SLAVE %2d is up" % rank)

    max_iters, c_rate, f_factor, pop_size, plen, pname = comm.bcast(None)

    rpack = rosetta_pack.RosettaPack(pname)

    pop = []

    def rand1bin():
        p1 = random.randint(0, pop_size - 1)
        p2 = random.randint(0, pop_size - 1)
        p3 = random.randint(0, pop_size - 1)

        while p1 == p2 or p2 == p3 or p1 == p3:
            p1 = random.randint(0, pop_size - 1)
            p2 = random.randint(0, pop_size - 1)
            p3 = random.randint(0, pop_size - 1)

        cutPoint = random.randint(0, plen)

        t_phi = []
        t_psi = []
        t_omega = []

        ind1 = pop[p1]
        ind2 = pop[p2]
        ind3 = pop[p3]

        for d in range(plen):
            if random.random() < c_rate or d == cutPoint:
                t_phi.append(ind1.phi[d] + (f_factor * (ind2.phi[d] - ind3.phi[d])))
                t_psi.append(ind1.psi[d] + (f_factor * (ind2.psi[d] - ind3.psi[d])))
                t_omega.append(ind1.omega[d] + (f_factor * (ind2.omega[d] - ind3.omega[d])))
            else:
                t_phi.append(ind1.phi[d])
                t_psi.append(ind1.psi[d])
                t_omega.append(ind1.omega[d])

            # if t_phi[d] != ind1.phi[d]:
                # print("%8.3f %8.3f" % (t_phi[d], ind1.phi[d]))

        trial = protein_data.ProteinData(rpack)
        trial.new_angles(t_phi, t_psi, t_omega)
        trial.fix_bounds()
        trial.eval()
        return trial

    # chunk_range = pop_size // (size - 1)

    # diff = 0
    # if rank == size - 1 and chunk_range * (size - 1) != pop_size:
    #     diff = pop_size - (chunk_range * (size - 1))

    # print('SLAVE %2d waking up' % rank)

    tag = 'WORK'
    comm.send(('BOOT', None), dest=0)

    while tag == 'WORK':
        tag, ind1, ind2, ind3, cutPoint = comm.recv()
        # print('SLAVE %2d got data' % rank)
        if tag == 'DIE':
            break

        t_phi = []
        t_psi = []
        t_omega = []

        for d in range(plen - 1):
            if random.random() < c_rate or d == cutPoint:
                t_phi.append(ind1.phi[d] + (f_factor * (ind2.phi[d] - ind3.phi[d])))
                t_psi.append(ind1.psi[d] + (f_factor * (ind2.psi[d] - ind3.psi[d])))
                t_omega.append(ind1.omega[d] + (f_factor * (ind2.omega[d] - ind3.omega[d])))
            else:
                t_phi.append(ind1.phi[d])
                t_psi.append(ind1.psi[d])
                t_omega.append(ind1.omega[d])

        trial = protein_data.ProteinData(rpack)
        trial.new_angles(t_phi, t_psi, t_omega)
        trial.fix_bounds()
        trial.eval()

        # print('SLAVE %2d sending data with score %8.2f' % (rank, trial.score))
        comm.send(('OK', trial), dest=0)
