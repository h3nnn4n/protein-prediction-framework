import config_loader
import sys
import de


def boot(conf_file):
    d = None

    cf = config_loader.ConfigLoader(conf_file)

    d = de.DE(pop_size=cf['pop_size'], max_iters=cf['max_iters'], pname=cf['pname'],
              f_factor=cf['f_factor'], c_rate=cf['c_rate'])

    cf.inject(d)

    d.reload_config()
    d.run()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        conf_file = sys.argv[1]
    else:
        conf_file = None

    boot(conf_file)
