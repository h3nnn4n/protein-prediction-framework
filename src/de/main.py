import config_loader
import comm
import sys
import de


if __name__ == '__main__':
    c = comm.Pigeon()
    d = None

    if len(sys.argv) > 1:
        conf_file = sys.argv[1]
    else:
        conf_file = None

    cf = config_loader.ConfigLoader(conf_file)

    d = de.DE(pop_size=cf['pop_size'], max_iters=cf['max_iters'], pname=cf['pname'],
              f_factor=cf['f_factor'], c_rate=cf['c_rate'], allatom=cf['allatom'])

    island = d
    island.set_coms(c)
    island.island_interval = 75

    island.run()
