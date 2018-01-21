import numpy as np
import boxplot
import hashlib
import sys
import re
import os


def hashfile(path, blocksize=65536):
    afile = open(path, 'rb')
    hasher = hashlib.md5()
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    afile.close()
    return hasher.hexdigest()


target = sys.argv[1]

os.chdir(target)

names = os.listdir()

hashs = set()

pnames = set()

# Organize the folders
if True:
    for name in names:
        if 'stats' in name:
            prefix = name[5:]

            stats_path = 'stats' + prefix
            parameters_path = 'parameters' + re.sub('.dat', '.yaml', prefix)

            base_name = prefix[2:6]

            if base_name not in pnames and not os.path.exists(base_name):
                os.mkdir(base_name)
                pnames.add(base_name)

            os.rename(stats_path, base_name + '/' + stats_path)
            os.rename(parameters_path, base_name + '/' + parameters_path)

            for n in os.listdir():
                if 'pdb' in n and prefix.split('.')[0] in n:
                    # print(prefix, n)
                    os.rename(n, base_name + '/' + n)

# Group similar files
if True:
    dnames = os.listdir()
    for dname in dnames:
        if len(dname) == 4:
            os.chdir(dname)

            names = os.listdir()
            for name in names:
                if 'stats' in name:

                    prefix = name[5:]

                    stats_path = 'stats' + prefix
                    parameters_path = 'parameters' + re.sub('.dat', '.yaml', prefix)

                    base_name = prefix[2:6]

                    h = hashfile(parameters_path)

                    p = h

                    # print(p, name)
                    if '____' not in name:
                        p = name.split('__')[2]
                        # print(p)

                    if h not in hashs and not os.path.exists(p):
                        os.mkdir(p)

                    hashs.add(h)

                    os.rename(stats_path, p + '/' + stats_path)
                    os.rename(parameters_path, p + '/' + parameters_path)

                    if not os.path.exists(p + '/pdb'):
                        os.mkdir(p + '/pdb')

                    for n in os.listdir():
                        if 'pdb' in n and prefix.split('.')[0] in n:
                            # print(prefix, n)
                            os.rename(n, p + '/pdb/' + n)

            os.chdir('..')

# Stats
if True:
    dnames = os.listdir()
    raws = {}
    for dname in dnames:
        if len(dname) == 4:
            if dname not in raws.keys():
                raws[dname] = []
            os.chdir(dname)

            alldata = {}

            names = os.listdir()
            for name in names:
                if not os.path.isdir(name):
                    continue

                os.chdir(name)
                cases = os.listdir()

                data = {}
                data['best'] = []
                data['mean'] = []
                # data['dist'] = []
                data['mdf'] = []
                data['rmsd'] = []

                for case in cases:
                    if 'stats' not in case:
                        continue

                    with open(case, 'rt') as f:
                        for l in f.readlines():
                            line = l.rstrip().lstrip()
                            tokens = re.sub("\s+", " ", line).split(' ')

                        if len(tokens) < 3:
                            # print(tokens)
                            continue

                        best = float(tokens[2])
                        mean = float(tokens[3])
                        mdf = float(tokens[4])
                        rmsd = float(tokens[6])
                        # dist = float(tokens[4])

                        data['best'].append(best)
                        data['mean'].append(mean)
                        # data['dist'].append(dist)
                        data['mdf'].append(mdf)
                        data['rmsd'].append(rmsd)

                    # print(case, data['rmsd'][0])

                alldata[name] = data
                raws[dname].append(name + '_raw.dat')
                with open('../' + name + '_raw.dat', 'wt') as f:
                    for i in range(len(data['best'])):
                        f.write('%f %f %f %f\n' % (data['best'][i], data['mean'][i], data['mdf'][i], data['rmsd'][i]))

                with open('../' + name + '_mean_std.dat', 'wt') as f:
                    f.write('%10.5f %10.5f %10.5f %10.5f\n' %
                            (np.mean(data['best']), np.mean(data['mean']), np.mean(data['mdf']), np.mean(data['rmsd'])))

                    f.write('%10.5f %10.5f %10.5f %10.5f\n' %
                            (np.std(data['best']), np.std(data['mean']), np.std(data['mdf']), np.std(data['rmsd'])))

                    f.write('%10.5f %10.5f %10.5f %10.5f\n' %
                            (np.median(data['best']), np.median(data['mean']), np.median(data['mdf']), np.median(data['rmsd'])))

                os.chdir('..')
            ps = ['best', 'mean', 'mdf', 'rmsd']

            with open('all_data.dat', 'wt') as f:
                for k, v in alldata.items():
                    for p in ps:
                        f.write('%25s %6s %f %f %f %f\n' % (k, p, min(v[p]), np.mean(v[p]), np.std(v[p]), np.median(v[p])))

            for p in ps:
                with open('all_%s.dat' % p, 'wt') as f:
                    for k, v in alldata.items():
                        f.write('%25s %6s %f %f %f %f\n' % (k, p, min(v[p]), np.mean(v[p]), np.std(v[p]), np.median(v[p])))

            os.chdir('..')

    for k, v in raws.items():
        os.chdir(k)
        boxplot.boxplot(k, v)
        os.chdir('..')
