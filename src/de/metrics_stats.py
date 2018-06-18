import pickle
import numpy as np


with open('data.pickle', 'rb') as f:
    p = pickle.load(f)

data, prots, tests, tests_prots = p

metrics = ['gdt_ts', 'gdt_ha', 'tm_score', 'maxsub']
modes = ['before', 'after']

general = {}

for mode in modes:
    general[mode] = {}
    for p in prots:
        general[mode][p] = {}
        for t in tests:
            general[mode][p][t] = {}
            for m in metrics:
                general[mode][p][t][m] = []

for prot in data:
    pdata = data[prot]

    for test in pdata:
        tdata = pdata[test]
        test_name = test[5:]

        for mode in tdata:
            runs = tdata[mode]

            for run in runs:
                general[mode][prot][test_name]['gdt_ts'].append(run['gdt_ts'][0])
                general[mode][prot][test_name]['gdt_ha'].append(run['gdt_ha'][0])
                general[mode][prot][test_name]['maxsub'].append(run['maxsub'])
                general[mode][prot][test_name]['tm_score'].append(run['tm_score'])

for mode in modes:
    for p in prots:
        for t in tests:
            for m in metrics:
                if t == 'sade_remc':
                    print(mode, p, m, np.mean(general[mode][p][t][m]))
        print()
    print()
