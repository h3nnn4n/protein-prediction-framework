import comm
import de


if __name__ == '__main__':
    c = comm.Pigeon()
    d = None

    pop_size = 100
    max_iters = 5000
    pname = '1zdd'
    pname = '1rop'
    pname = '1crn'
    # pname = '1plw'

    c_rate = 1.0
    f_factor = 0.5

    allatom = False

    # if c.rank % 2 == 0:
        # d = de.DE(pop_size=pop_size, max_iters=max_iters, pname=pname, allatom=True)
    # else:
        # d = de.DE(pop_size=pop_size, max_iters=max_iters, pname=pname)

    # d = de.DE(pop_size=pop_size, max_iters=max_iters, pname=pname, allatom=True)
    d = de.DE(pop_size=pop_size, max_iters=max_iters, pname=pname, f_factor=f_factor, c_rate=c_rate, allatom=allatom)

    island = d
    island.set_coms(c)
    island.island_interval = 75
    # island.island_interval = 2

    # if island.comm.rank == 1:
        # island.coil_only = True

    # if island.comm.rank == 2:
        # island.c_rate = 0.3

    island.run()
