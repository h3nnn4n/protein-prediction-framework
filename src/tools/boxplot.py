import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats
import scipy
import numpy as np
import sys
import os
import re


def boxplot(target, names):
    alldata = {}
    for name in names:
        if os.path.isdir(name) or 'raw' not in name or target not in name:
            continue

        data = {}
        data['best'] = []
        data['mean'] = []
        data['mdf'] = []
        data['best_fxn'] = []
        data['rmsd'] = []
        data['rmsd_fxn'] = []

        t = name

        with open(t, 'rt') as f:
            for l in f.readlines():
                line = l.rstrip().lstrip()
                tokens = re.sub("\s+", " ", line).split(' ')

                if len(tokens) < 2:
                    continue

                best = float(tokens[0])
                best_fxn = float(tokens[1])
                mean = float(tokens[2])
                mdf = float(tokens[3])
                rmsd = float(tokens[4])
                rmsd_fxn = float(tokens[5])

                data['best'].append(best)
                data['best_fxn'].append(best_fxn)
                data['mean'].append(mean)
                data['mdf'].append(mdf)
                data['rmsd'].append(rmsd)
                data['rmsd_fxn'].append(rmsd_fxn)

        alldata[name] = data

    x = []
    names = []
    keys = []
    for k, v in alldata.items():
        x.append(v['best_fxn'])
        # x.append(v['rmsd'])
        names.append(k.split('.')[0][:-4])
        keys.append((k, k.split('.')[0][:-4]))

    names.sort()

    # print(alldata)

    for mode in ['rmsd_fxn', 'best_fxn']:
        for a, i in enumerate(keys):
            for b, j in enumerate(keys[a + 1:]):
                p1 = alldata[i[0]][mode]
                p2 = alldata[j[0]][mode]

                ml = min(len(p1), len(p2))

                w = scipy.stats.wilcoxon(p1[:ml], p2[:ml])[1]

                if w < 0.05:
                    a, b = np.mean(p1[:ml]), np.mean(p2[:ml])
                    if a < b:
                        print("%8s %35s %35s %8.5f %8.3f %8.3f %8d" % (mode, i[1], j[1], w, a, b, ml))
                    else:
                        print("%8s %35s %35s %8.5f %8.3f %8.3f %8d" % (mode, j[1], i[1], w, b, a, ml))

    try:
        ############

        names = sorted([k.split('.')[0][:-4] for k in alldata.keys()])
        values_dict = {k.split('.')[0][:-4]: v['best_fxn'] for k, v in alldata.items()}
        x = [values_dict[name] for name in names]

        fig, ax = plt.subplots()

        ax.set_title(target.upper() + ' Boxplot for best Energy')
        # ax.set_xlabel('xlabel')
        ax.set_ylabel('Energy')

        ax.boxplot(x)

        xtickNames = plt.setp(ax, xticklabels=np.repeat(names, 1))
        plt.setp(xtickNames, rotation=90, fontsize=10)
        plt.tight_layout()

        plt.savefig(target + '_energy_boxplot.png')
        plt.close()

        ############

        names = sorted([k.split('.')[0][:-4] for k in alldata.keys()])
        values_dict = {k.split('.')[0][:-4]: v['rmsd_fxn'] for k, v in alldata.items()}
        x = [values_dict[name] for name in names]

        fig, ax = plt.subplots()

        ax.set_title(target.upper() + ' Boxplot for best RMSD')
        # ax.set_xlabel('xlabel')
        ax.set_ylabel('RMSD')

        ax.boxplot(x)

        xtickNames = plt.setp(ax, xticklabels=np.repeat(names, 1))
        plt.setp(xtickNames, rotation=90, fontsize=10)
        plt.tight_layout()

        plt.savefig(target + '_rmsd_boxplot.png')
        plt.close()
    except Exception:
        print('Could not render plots')
        pass


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('insuficient parameters')

    target = sys.argv[1]
    names = sys.argv[2:]

    boxplot(target, names)
