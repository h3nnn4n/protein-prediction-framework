import sys
import os
import os.path
import pickle

sys.path.append('../../external')

import tmscore.src.TMscore as tmscore


tmscore = tmscore.TMscore(path='../../external/tmscore/src/TMscore')

path = {}
path['1zdd'] = '/home/h3nnn4n/progs/de_supimpa/protein_data/1zdd/1zdd.pdb'
path['1crn'] = '/home/h3nnn4n/progs/de_supimpa/protein_data/1crn/1crn.pdb'
path['1enh'] = '/home/h3nnn4n/progs/de_supimpa/protein_data/1enh/1enh.pdb'
path['1ail'] = '/home/h3nnn4n/progs/de_supimpa/protein_data/1ail/1ail.pdb'

data = {}
prots = set()
tests = set()
tests_prots = set()


def find_energy(name):
    prefix = name

    ind = name.find('__') + 2
    prefix = name[ind:-4]

    cwd = os.getcwd()
    os.chdir('..')

    names = os.listdir()

    for name in names:
        if prefix in name:
            with open(name) as f:
                for line in f.readlines():
                    if 'scorefxn' in line:
                        values = line.split(':')
                        energy = float(values[1].strip())

    os.chdir(cwd)

    return energy


def do_stuff():
    dnames = os.listdir()
    for dname in dnames:
        if len(dname) == 4:
            # print(dname)
            prots.add(dname)
            if dname not in data.keys():
                data[dname] = {}
            os.chdir(dname)

            names = os.listdir()
            for name in names:
                if not os.path.isdir(name):
                    continue

                tests_prots.add(name)
                tests.add(name[5:])

                if name not in data[dname].keys():
                    data[dname][name] = {}

                # print(' ', name)
                os.chdir(name)
                os.chdir('pdb')
                pdbs = os.listdir()

                for pdb in pdbs:
                    tokens = pdb.split('__')
                    mode = 'before' if 'repacked' in tokens[0] else 'after'
                    if mode not in data[dname][name].keys():
                        data[dname][name][mode] = []

                    # print('  ', pdb)
                    tmscore(path[dname], pdb)
                    scores = tmscore.get_all()
                    scores['scorefxn'] = find_energy(pdb)
                    data[dname][name][mode].append(scores)

                os.chdir('../..')

            os.chdir('..')


if __name__ == '__main__':
    do_stuff()
    ninja = (data, prots, tests, tests_prots)
    with open('data.pickle', 'wb') as f:
        pickle.dump(ninja, f)
