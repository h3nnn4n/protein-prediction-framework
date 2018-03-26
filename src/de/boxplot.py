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
        data['rmsd'] = []

        t = name

        with open(t, 'rt') as f:
            for l in f.readlines():
                line = l.rstrip().lstrip()
                tokens = re.sub("\s+", " ", line).split(' ')

                if len(tokens) < 2:
                    continue

                best = float(tokens[0])
                mean = float(tokens[1])
                mdf = float(tokens[2])
                rmsd = float(tokens[3])

                data['best'].append(best)
                data['mean'].append(mean)
                data['mdf'].append(mdf)
                data['rmsd'].append(rmsd)

        alldata[name] = data

    x = []
    names = []
    keys = []
    for k, v in alldata.items():
        x.append(v['best'])
        # x.append(v['rmsd'])
        names.append(k.split('.')[0][:-4])
        keys.append((k, k.split('.')[0][:-4]))

    names.sort()

    # print(alldata)

    for mode in ['rmsd', 'best']:
        for a, i in enumerate(keys):
            for b, j in enumerate(keys[a + 1:]):
                p1 = alldata[i[0]][mode]
                p2 = alldata[j[0]][mode]

                ml = min(len(p1), len(p2))

                w = scipy.stats.wilcoxon(p1[:ml], p2[:ml])[1]

                if w < 0.05:
                    print("%8s %26s %26s %8.5f %8.5f %8.5f %8d" % (mode, i[1], j[1], w, np.mean(p1[:ml]), np.mean(p2[:ml]), ml))

    fig, ax = plt.subplots()

    # ax.set_title('1ZDD Boxplot for best Energy')
    # ax.set_xlabel('xlabel')
    ax.set_ylabel('Energy')

    ax.boxplot(x)

    xtickNames = plt.setp(ax, xticklabels=np.repeat(names, 1))
    plt.setp(xtickNames, rotation=90, fontsize=10)
    plt.tight_layout()

    plt.savefig(target + '_energy_boxplot.png')

    x = []
    names = []
    for k, v in alldata.items():
        x.append(v['rmsd'])
        # x.append(v['rmsd'])
        names.append(k.split('.')[0][:-4])
    names.sort()

    fig, ax = plt.subplots()

    # ax.set_title('1ZDD Boxplot for best RMSD')
    # ax.set_xlabel('xlabel')
    ax.set_ylabel('RMSD')

    ax.boxplot(x)

    xtickNames = plt.setp(ax, xticklabels=np.repeat(names, 1))
    plt.setp(xtickNames, rotation=90, fontsize=10)
    plt.tight_layout()

    plt.savefig(target + '_rmsd_boxplot.png')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('insuficient parameters')

    target = sys.argv[1]
    names = sys.argv[2:]

    boxplot(target, names)
