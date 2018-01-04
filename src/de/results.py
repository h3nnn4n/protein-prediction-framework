class Results:
    def __init__(self, stats, config):
        self.stats = stats
        self.config = config

        self.stats_data = self.load(stats)
        self.config_data = self.load(config)

    def load(self, name):
        data = []
        with open(name) as f:
            for l in f.readlines():
                line = l.rstrip().lstrip()
                data.append(line)

        return data

    def store(self):
        def dump(name, data):
            with open(name, 'wt') as f:
                for l in data:
                    f.write(l + '\n')

        dump(self.stats, self.stats_data)
        dump(self.config, self.config_data)
