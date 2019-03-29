import yaml
import os


class ConfigLoader:
    def __init__(self, conf_file=None):
        self.conf_file = conf_file
        self.set_default_options()

        try:
            self.inject_from_file()
        except Exception as e:
            print('Using default values because: %s\nWhile opening %s' % (e, conf_file))
            if conf_file is None:
                print('conf_file was not set')
            else:
                if os.path.isfile(conf_file):
                    print('file exists')
                else:
                    print('file does not exist')
                    for p in conf_file:
                        print(p)
            self.p_values = self.defaults

    def __getitem__(self, key):
        return self.p_values[self.parameters.index(key)]

    def inject(self, de):
        for k, v in zip(self.parameters, self.p_values):
            if k not in de.__dict__.keys():
                raise KeyError('WARN: key %s was not found!' % k)

            de.__dict__[k] = v

    def inject_from_file(self):
        config = self.get_config_file()

        for n, (k, _) in enumerate(zip(self.parameters, self.defaults)):
            if k in config.keys():
                self.p_values.append(config[k])
            else:
                self.p_values.append(self.defaults[n])

    def get_config_file(self):
        with open(self.conf_file, 'r') as f:
            config = yaml.load(f)

        return config

    def set_default_options(self):
        conf_file = self.conf_file

        self.options = {}

        self.options['pname'] = '1crn'
        self.options['pop_size'] = 10
        self.options['max_iters'] = 10
        self.options['max_evals'] = 500
        self.options['stop_condition'] = 'evals'
        self.options['c_rate'] = 1.0
        self.options['f_factor'] = 0.5
        self.options['mode'] = 'normal'

        self.options['stage0_init'] = False
        self.options['stage2_interval'] = -1
        self.options['stage2_all_interval'] = -1
        self.options['partial_reset'] = -1

        self.options['log_interval'] = 10

        self.options['do_lsh'] = False
        self.options['n_hashes'] = 10
        self.options['n_buckets'] = 50
        self.options['update_interval'] = 20
        self.options['change_interval'] = 100

        self.options['reset_d_trigger'] = 0.0
        self.options['reset_d_percent'] = 0.0
        self.options['reset_rmsd_trigger'] = 0.0
        self.options['reset_rmsd_percent'] = 0.0

        self.options['forced_insertion'] = False
        self.options['forced_insertion_chance'] = 0.01
        self.options['forced_insertion_mode'] = 'frag9'

        self.options['cname'] = 'none' if conf_file is None else conf_file.split('.')[0]
        self.options['sade_run'] = False
        self.options['sade_lp'] = 50
        self.options['sade_reinit_interval'] = 1000
        self.options['sade_selection'] = 'roulette'

        self.options['enable_remc'] = True

        self.options['do_crowding'] = False
        self.options['do_rmsd_crowding'] = False
        self.options['crowding_factor'] = 3

        self.options['ops'] = ['rand1bin_global']

        self.options['do_clearing'] = False
        self.options['clearing_interval'] = 10
        self.options['clearing_size'] = 1

        self.options['energy_function'] = 'score3'
        self.options['energy_options'] = ['score0_200', 'score1_400', 'score2_1400', 'score3_2000', 'score5_1000']

        self.options['extended_diversity_measurements'] = True

        self.parameters = []
        self.defaults = []
        self.p_values = []

        for k in self.options.keys():
            self.parameters.append(k)
            self.defaults.append(self.options[k])
