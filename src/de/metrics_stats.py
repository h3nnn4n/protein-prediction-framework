import pickle
import numpy as np
import pandas as pd
import matplotlib as mpl
import pandas.tools.plotting

from sklearn import preprocessing

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


with open('data.pickle', 'rb') as f:
    p = pickle.load(f)

data, prots, tests, tests_prots = p

metrics = ['scorefxn', 'gdt_ts', 'gdt_ha', 'tm_score', 'maxsub', 'rmsd']
#metrics = ['gdt_ts', 'gdt_ha', 'tm_score', 'maxsub', 'rmsd']
#metrics = ['gdt_ts', 'gdt_ha', 'tm_score', 'maxsub']
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

values = {}

#target_protein = '1zdd'
#target_protein = '1crn'
target_protein = '1enh'
target_protein = '1ail'

for mode in modes:
    if mode != 'before':
        continue
    for p in prots:
        if target_protein != p:
            continue
        for t in tests:
            values[t] = []
            print()
            print(t)
            #for _ in metrics:
                #for _ in range(10):
            #for _ in range(10):
                    #labels.append(t)
                #labels.append(t)
            labels.append(t)

            for m in metrics:
                v = (general[mode][p][t][m])

                #v = np.mean(v)
                #if m in ['scorefxn', 'rmsd']:
                    #v = min(v)
                #else:
                    #v = max(v)

                #print(v)
                values[t].append(v)

                # print(t)
                #if t != 'sade_remc':
                    #continue

                #print(v)
                #continue

                m_ = True

                if p == target_protein and m == 'tm_score':
                    if m_:
                        tm.append(max(general[mode][p][t][m]))
                    else:
                        tm.append(general[mode][p][t][m])
                if p == target_protein and m == 'rmsd':
                    #rmsd.append(general[mode][p][t][m])
                    #print(general[mode][p][t][m])
                    if m_:
                        rmsd.append(min(general[mode][p][t][m]))
                    else:
                        rmsd.append(general[mode][p][t][m])
                if p == target_protein and m == 'gdt_ha':
                    if m_:
                        gdt_ha.append(max(general[mode][p][t][m]))
                    else:
                        gdt_ha.append(general[mode][p][t][m])
                if p == target_protein and m == 'gdt_ts':
                    if m_:
                        gdt_ts.append(max(general[mode][p][t][m]))
                    else:
                        gdt_ts.append(general[mode][p][t][m])
                if p == target_protein and m == 'maxsub':
                    if m_:
                        maxsub.append(max(general[mode][p][t][m]))
                    else:
                        maxsub.append(general[mode][p][t][m])
                if p == target_protein and m == 'scorefxn':
                    if m_:
                        scorefxn.append(min(general[mode][p][t][m]))
                    else:
                        scorefxn.append(general[mode][p][t][m])
        # print()
    # print()

#n_groups = 8
#n_groups = len(metrics)

#fig, ax = plt.subplots()

#index = np.arange(n_groups)
#bar_width = 0.1

#opacity = 0.4
#error_config = {'ecolor': '0.3'}

#bars = []

#print(n_groups, len(index), len(values[list(tests)[0]]))

#print(values[list(tests)[0]])
#plt.bar(index, values[list(tests)[0]])

#for k, t in enumerate(tests):
#for m in metrics:
    #a = ax.bar(index + k * bar_width, values[t], bar_width,
               #alpha=opacity)

    #bars.append(a)
    #continue

#plt.savefig(target_protein + '_bars.png')

#for k, v in values.items():
    #print(k, v)
    #print()

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

flatten = lambda l: [item for sublist in l for item in sublist]
#flatten = lambda i: i

min_max_scaler = preprocessing.MinMaxScaler()

def scale(a):
    #b = np.asarray(flatten(a)).reshape(-1, 1)
    b = np.asarray((a)).reshape(-1, 1)
    c = min_max_scaler.fit_transform(b)
    c = flatten(c.reshape(-1, 1))

    return c

#kappa = {
            #'gdt_ts': scale(gdt_ts),
            #'gdt_ha': scale(gdt_ha),
            #'tm_score': scale(tm),
            #'maxsub': scale(maxsub),
            #'rmsd': scale(rmsd),
            #'scorefxn': scale(scorefxn),
            #'labels': labels
        #}

#kappa = {
            #'gdt_ts': flatten(gdt_ts),
            #'gdt_ha': flatten(gdt_ha),
            #'tm_score': flatten(tm),
            #'maxsub': flatten(maxsub),
            #'rmsd': flatten(rmsd),
            #'scorefxn': flatten(scorefxn),
            #'labels': labels
        #}

def invert(a):
    return list(map(lambda x: (1.0 - x) * 0.15 + 0.25, a))

kappa = {
            'gdt_ts': (gdt_ts),
            'gdt_ha': (gdt_ha),
            'tm_score': (tm),
            'maxsub': (maxsub),
            #'rmsd': invert(scale(rmsd)),
            #'scorefxn': invert(scale(scorefxn)),
            'labels': labels
        }

#print(kappa)

#for k, v in kappa.items():
    #print(k, len(v))
    #print(k, (v))

df = pd.DataFrame.from_dict(kappa, orient='columns')

#print(df.head())

#grid = sns.PairGrid(data=df)
##grid = grid.map_upper(plt.scatter)
#grid = grid.map_upper(sns.regplot)
#grid = grid.map_upper(corr)
#grid = grid.map_diag(plt.hist)
##grid = grid.map_diag(sns.kdeplot)
#grid = grid.map_lower(sns.kdeplot)

#plt.savefig(target_protein + '_pairplot.png')

sns.color_palette("PuBuGn_d")

f = plt.figure(figsize=(8, 5))
pandas.tools.plotting.parallel_coordinates(df, 'labels')
lgd = plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

plt.title(target_protein + ' parallel plot', color='black')

plt.tight_layout()

plt.savefig(target_protein + '_parallel.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
