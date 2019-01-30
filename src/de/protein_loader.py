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

        original = os.path.dirname(os.path.realpath(__file__))
        self.original = os.getcwd()

        base = self.protein_data_path

        dest = base + '/' + self.name
        if os.path.isdir(dest):
            os.chdir(dest)
            base_file = name + '.pdb'
            native_path = os.path.join(base, name, base_file)
            if os.path.isfile(native_path):
                self.native_path = native_path
            else:
                raise FileNotFoundError('Could not find %s' % native_path)

            if os.path.isfile(name + '.fasta'):
                with open(name + '.fasta', 'rt') as f:
                    self.target = ''
                    for line in f.readlines()[1:]:
                        l = line.rstrip()
                        self.target += l

                self.target = self.target.split('>')[0]

            if os.path.isfile(name + '.psipred.ss2'):
                with open(name + '.psipred.ss2', 'rt') as f:
                    self.ss_pred = ''
                    for line in f.readlines()[1:]:
                        tokens = re.sub(r"\s+", " ", line.rstrip().lstrip()).split(' ')
                        if len(tokens) < 3:
                            continue

                        self.ss_pred += tokens[2]

            f3 = base + '/' + name + '/output/' + name + '.200.3mers'
            if os.path.isfile(f3):
                self.fragset3_path = f3
            else:
                raise FileNotFoundError('Could not find %s' % f3)

            f9 = base + '/' + name + '/output/' + name + '.200.9mers'
            if os.path.isfile(f9):
                self.fragset9_path = f9
            else:
                raise FileNotFoundError('Could not find %s' % f9)
        else:
            raise FileNotFoundError('Could not find base folder: %s' % dest)

        os.chdir(original)

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
