import pickle
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


with open('data.pickle', 'rb') as f:
    p = pickle.load(f)

data, prots, tests, tests_prots = p

metrics = ['gdt_ts', 'gdt_ha', 'tm_score', 'maxsub', 'rmsd']
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
                general[mode][prot][test_name]['rmsd'].append(run['rmsd'])

rmsd = []
maxsub = []
gdt_ts = []
gdt_ha = []
tm = []
labels = []

for mode in modes:
    if mode != 'before':
        continue
    for p in prots:
        for t in tests:
            for m in metrics:
                # if p == '1zdd' and m == 'gdt_ts':
                # if p == '1crn' and m == 'gdt_ts':
                # if p == '1crn' and m == 'maxsub':
                if p == '1zdd' and m == 'tm_score':
                    tm.append(general[mode][p][t][m])
                    labels.append(t)
                if p == '1zdd' and m == 'rmsd':
                    rmsd.append(general[mode][p][t][m])
                if p == '1zdd' and m == 'gdt_ha':
                    gdt_ha.append(general[mode][p][t][m])
                if p == '1zdd' and m == 'gdt_ts':
                    gdt_ts.append(general[mode][p][t][m])
                if p == '1zdd' and m == 'maxsub':
                    maxsub.append(general[mode][p][t][m])
                # if t == 'sade_remc':
                    # print("%10s %4s %10s %8.4f" % (mode, p, m, np.mean(general[mode][p][t][m])))
                    # if m == 'rmsd':
                        # print("%35s %4s %10s %8.4f" % (t, p, m, min(general[mode][p][t][m])))
                    # else:
                        # print("%35s %4s %10s %8.4f" % (t, p, m, max(general[mode][p][t][m])))
        # print()
    # print()
labels=list('ABCDEFGH')
fs = 10
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 12), sharey=False)

axes[0, 0].boxplot(rmsd, labels=labels)
axes[0, 0].set_title('RMSD', fontsize=fs)

# axes[0, 0].boxplot(gdt_ha, labels=labels)
# axes[0, 0].set_title('GDT-HA', fontsize=fs)

axes[0, 1].boxplot(gdt_ts, labels=labels)
axes[0, 1].set_title('GDT-TS', fontsize=fs)

axes[1, 1].boxplot(tm, labels=labels)
axes[1, 1].set_title('TM-Score', fontsize=fs)

axes[1, 0].boxplot(maxsub, labels=labels)
axes[1, 0].set_title('MaxSub', fontsize=fs)

fig.subplots_adjust(hspace=0.4)


plt.savefig('1zdd_rmsd.png')
