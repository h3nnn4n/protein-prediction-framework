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
            repack_path = 'repack' + prefix
            ops_path = 'ops' + prefix
            parameters_path = 'parameters' + re.sub('.dat', '.yaml', prefix)

            base_name = prefix[2:6]

            if base_name not in pnames and not os.path.exists(base_name):
                os.mkdir(base_name)
                pnames.add(base_name)

            os.rename(stats_path, base_name + '/' + stats_path)
            os.rename(parameters_path, base_name + '/' + parameters_path)

            if os.path.exists(ops_path):
                os.rename(ops_path, base_name + '/' + ops_path)

            if os.path.exists(repack_path):
                os.rename(repack_path, base_name + '/' + repack_path)

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
                    repack_path = 'repack' + prefix
                    ops_path = 'ops' + prefix
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

                    if os.path.exists(ops_path):
                        os.rename(ops_path, p + '/' + ops_path)

                    if os.path.exists(repack_path):
                        os.rename(repack_path, p + '/' + repack_path)

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
                data['best_fxn'] = []
                # data['dist'] = []
                data['mdf'] = []
                data['rmsd'] = []
                data['rmsd_fxn'] = []

                for case in cases:
                    if 'stats' in case:
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
                    elif 'repack' in case:
                        with open(case, 'rt') as f:
                            for line in f.readlines():
                                tokens = re.sub("\s+", " ", line.strip()).split(' ')

                                if 'scorefxn' in tokens[0]:
                                    score = float(tokens[1])
                                    # print('scorefxn', d)

                                if 'rmsd_after' in tokens[0]:
                                    rmsd = float(tokens[1])
                                    # print('rmsd_after', d)

                            data['rmsd_fxn'].append(rmsd)
                            data['best_fxn'].append(score)

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

            ps = ['best_fxn', 'best', 'mean', 'mdf', 'rmsd', 'rmsd_fxn']

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
