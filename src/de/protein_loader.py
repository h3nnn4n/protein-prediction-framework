import os
import re


class ProteinLoader:
    def __init__(self, name=None):
        if name is not None:
            self.name = name

        self.fragset3_path = None
        self.fragset9_path = None

        self.target = None
        self.ss_pred = None

        self.native_path = None

        self.protein_data_path = '../../protein_data'

    def load(self, name=None):
        if name is not None:
            self.name = name

        base = self.protein_data_path

        dest = base + '/' + self.name
        if os.path.isdir(dest):
            os.chdir(dest)
            if os.path.isfile(name + '.pdb'):
                self.native_path = name + '.pdb'

            if os.path.isfile(name + '.fasta'):
                with open(name + '.fasta', 'rt') as f:
                    self.target = ''
                    for line in f.readlines()[1:]:
                        l = line.rstrip()
                        self.target += l

            if os.path.isfile(name + '.psipred.ss2'):
                with open(name + '.psipred.ss2', 'rt') as f:
                    self.ss_pred = ''
                    for line in f.readlines()[1:]:
                        tokens = re.sub("\s+", " ", line.rstrip().lstrip()).split(' ')
                        if len(tokens) < 3:
                            continue

                        self.ss_pred += tokens[2]

            f3 = base + '/' + name + '/output/' + name + '.200.3mers'
            if os.path.isfile(f3):
                self.fragset3_path = f3

            f9 = base + '/' + name + '/output/' + name + '.200.9mers'
            if os.path.isfile(f9):
                self.fragset9_path = f9

    def show(self):
        print(self.name)
        print(self.fragset3_path)
        print(self.fragset9_path)
        print(self.target)
        print(self.ss_pred)
        print(self.native_path)

    def get_data(self):
        return (self.name, self.target, self.ss_pred, self.native_path, self.fragset3_path, self.fragset9_path)


if __name__ == '__main__':
    x = ProteinLoader()
    x.load('1rop')
    x.show()
