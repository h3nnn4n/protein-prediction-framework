import yaml


class ConfigLoader:
    def __init__(self, conf_file):
        self.parameters = ['pname', 'pop_size', 'max_iters', 'c_rate', 'f_factor', 'allatom', 'coil_only', 'stage0_init', 'stage0_init',
                           'stage2_interval', 'stage2_all_interval', 'partial_reset', 'log_interval', 'island_interval', 'do_lhs',
                           'n_hashes', 'update_interval', 'change_interval']
        self.defaults = ['1crn', 100, 50, 1.0, 0.5, False, False, False, False, -1, -1, -1, 10, 100, False, 10, 20, 100]
        self.p_values = []

        try:
            with open(conf_file, 'r') as f:
                try:
                    config = yaml.load(f)

                    for n, (k, v) in enumerate(zip(self.parameters, self.defaults)):
                        if k in config.keys():
                            self.p_values.append(config[k])
                        else:
                            self.p_values.append(self.defaults[n])
                except yaml.YAMLError as ee:
                    print(ee)
        except Exception as e:
            # print(e)
            self.p_values = self.defaults

        # print(self.p_values)

    def __getitem__(self, key):
        return self.p_values[self.parameters.index(key)]

    def inject(self, de):
        for k, v in zip(self.parameters, self.p_values):
            de.__dict__[k] = v