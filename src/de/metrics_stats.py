import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl

mpl.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns


def corr(x, y, **kwargs):

    # Calculate the value
    coef = np.corrcoef(x, y)[0][1]
    # Make the label
    label = r'$\rho$ = ' + str(round(coef, 2))

    # Add the label to the plot
    ax = plt.gca()
    ax.annotate(label, xy=(0.2, 0.95), size=20, xycoords=ax.transAxes)

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


with open('data.pickle', 'rb') as f:
    p = pickle.load(f)

data, prots, tests, tests_prots = p

metrics = ['scorefxn', 'gdt_ts', 'gdt_ha', 'tm_score', 'maxsub', 'rmsd']
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
                general[mode][prot][test_name]['scorefxn'].append(run['scorefxn'])

tm = []
rmsd = []
maxsub = []
gdt_ts = []
gdt_ha = []
scorefxn = []
labels = []

target_protein = '1zdd'

for mode in modes:
    if mode != 'before':
        continue
    for p in prots:
        # print()
        for t in tests:
            # print(t)
            for m in metrics:
                # print(t)
                if t != 'sade_remc':
                    continue

                if p == target_protein and m == 'tm_score':
                    tm.append(general[mode][p][t][m])
                    labels.append(t)
                if p == target_protein and m == 'rmsd':
                    rmsd.append(general[mode][p][t][m])
                if p == target_protein and m == 'gdt_ha':
                    gdt_ha.append(general[mode][p][t][m])
                if p == target_protein and m == 'gdt_ts':
                    gdt_ts.append(general[mode][p][t][m])
                if p == target_protein and m == 'maxsub':
                    maxsub.append(general[mode][p][t][m])
                if p == target_protein and m == 'scorefxn':
                    scorefxn.append(general[mode][p][t][m])
        # print()
    # print()

# labels = list('ABCDEFGH')
# fs = 10
# fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16, 10), sharey=False)

# axes[0, 0].boxplot(scorefxn, labels=labels)
# axes[0, 0].set_title('scorefxn', fontsize=fs)

# axes[0, 1].boxplot(rmsd, labels=labels)
# axes[0, 1].set_title('RMSD', fontsize=fs)

# axes[0, 2].boxplot(tm, labels=labels)
# axes[0, 2].set_title('TM-Score', fontsize=fs)

# axes[1, 0].boxplot(maxsub, labels=labels)
# axes[1, 0].set_title('MaxSub', fontsize=fs)

# axes[1, 1].boxplot(gdt_ha, labels=labels)
# axes[1, 1].set_title('GDT-HA', fontsize=fs)

# axes[1, 2].boxplot(gdt_ts, labels=labels)
# axes[1, 2].set_title('GDT-TS', fontsize=fs)

# fig.subplots_adjust(hspace=0.4)

# plt.savefig(target_protein + '.png')

# plt.close()


kappa = {
            'scorefxn': scorefxn[0],
            'gdt_ts': gdt_ts[0],
            'gdt_ha': gdt_ha[0],
            'tm_score': tm[0],
            'maxsub': maxsub[0],
            'rmsd': rmsd[0]
        }

df = pd.DataFrame.from_dict(kappa, orient='columns')

grid = sns.PairGrid(data=df)
#grid = grid.map_upper(plt.scatter)
grid = grid.map_upper(sns.regplot)
grid = grid.map_upper(corr)
grid = grid.map_diag(plt.hist)
#grid = grid.map_diag(sns.kdeplot)
grid = grid.map_lower(sns.kdeplot)

plt.savefig(target_protein + '_pairplot.png')
