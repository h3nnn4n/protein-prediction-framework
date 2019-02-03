import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats
import pickle
import scipy
import numpy as np
import sys
import os
import re


def data_dump(target, names):
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

    with open('%s.pickle' % target, 'wb') as f:
        pickle.dump(alldata, f)
