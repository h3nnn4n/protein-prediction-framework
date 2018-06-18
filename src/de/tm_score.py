import sys
import os
import os.path

sys.path.append('../external')

import tmscore.src.TMscore as tmscore


tmscore = tmscore.TMscore(path='../external/tmscore/src/TMscore')

path = {}
path['1zdd'] = '/home/h3nnn4n/progs/de_supimpa/protein_data/1zdd/1zdd.pdb'
path['1crn'] = '/home/h3nnn4n/progs/de_supimpa/protein_data/1crn/1crn.pdb'
path['1enh'] = '/home/h3nnn4n/progs/de_supimpa/protein_data/1enh/1enh.pdb'
path['1ail'] = '/home/h3nnn4n/progs/de_supimpa/protein_data/1ail/1ail.pdb'

data = {}
prots = set()
tests = set()
tests_prots = set()


def a():
    dnames = os.listdir()
    for dname in dnames:
        if len(dname) == 4:
            print(dname)
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

                print(' ', name)
                os.chdir(name)
                os.chdir('pdb')
                pdbs = os.listdir()

                for pdb in pdbs:
                    tokens = pdb.split('__')
                    mode = 'before' if 'repacked' in tokens[0] else 'after'
                    if mode not in data[dname][name].keys():
                        data[dname][name][mode] = []

                    # print('  ', mode, ' ', pdb.split('__'))
                    print('  ', pdb)
                    tmscore(path[dname], pdb)
                    scores = tmscore.get_all()
                    data[dname][name][mode].append(scores)

                os.chdir('../..')

            os.chdir('..')


a()
