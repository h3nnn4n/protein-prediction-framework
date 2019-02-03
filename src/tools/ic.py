import scipy.stats
import numpy as np
import sys
import os
import re


final_data = {}


def ic(target, names):
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
                tokens = re.sub(r"\s+", " ", line).split(' ')

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

    # import code; code.interact(local=dict(globals(), **locals()))

    names.sort()

    for mode in ['rmsd_fxn', 'best_fxn']:
        if mode not in final_data.keys():
            final_data[mode] = {}

        x = final_data[mode]

        for i in keys:
            k = i[1][5:]
            if k not in x.keys():
                x[k] = []

            p1 = alldata[i[0]][mode]
            x[k].extend(p1)


def wrap():
    dnames = os.listdir()
    raws = {}
    for dname in dnames:
        if len(dname) == 4:
            if dname not in raws.keys():
                raws[dname] = []
            os.chdir(dname)

            names = os.listdir()
            for name in names:
                if not os.path.isdir(name):
                    continue

                raws[dname].append(name + '_raw.dat')

            os.chdir('..')

    for k, v in raws.items():
        os.chdir(k)
        ic(k, v)
        os.chdir('..')


if __name__ == '__main__':
    wrap()

    for k in final_data.keys():
        print()
        for l in final_data[k].keys():
            y = final_data[k][l]
            w = scipy.stats.t.interval(0.95, len(y) - 1, loc=np.mean(y), scale=scipy.stats.sem(y))
            print("%53s %8.4f %8.4f %8.4f" % (l, w[0], w[1], w[1] - w[0]))
