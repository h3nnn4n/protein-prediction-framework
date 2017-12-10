import matplotlib.pyplot as plt
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
    for k, v in alldata.items():
        x.append(v['best'])
        # x.append(v['rmsd'])
        names.append(k.split('.')[0][:-4])

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