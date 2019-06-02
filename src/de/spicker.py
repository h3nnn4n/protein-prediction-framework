import subprocess
import string
import random
import uuid
import os
import re


from shutil import rmtree


class Spicker:
    def __init__(self, de=None):
        self.de = de
        self.pop = self.de.pop
        self.run_spicker = self.de.run_spicker

        self.reset()

    def reset(self):
        self.working_folder = '%s__%s' % ('spicker', uuid.uuid4())

        self.tra1_filenames = []
        self.psize = None

        self.cluster_map = {}
        self.centroid_index_in_pop = []

    def run(self):
        if not self.run_spicker:
            return

        self.create_working_folder()
        self.init_spicker_files()
        self.remove_working_folder()

    def init_spicker_files(self):
        print('[SPICKER] Initializing')
        self.init_tra1_files()
        self.init_tra1in_file()
        self.init_rmsinp_file()
        self.init_seq_file()
        self.init_ca_file()

        print('[SPICKER] Running')
        self.call_spicker()

        print('[SPICKER] Parsing results')
        self.parse_results()

    def init_tra1_files(self):
        for index, individual in enumerate(self.pop):
            self.dump_individual_to_tra1(individual, index)

    def dump_individual_to_tra1(self, individual, index):
        filename = '%s.tra1' % random_string()
        out = os.path.join(self.working_folder, filename)
        pose = individual.pose
        score = individual.score

        self.cluster_map[index] = filename
        self.cluster_map[filename] = index

        data = []

        for i in range(pose.total_residue()):
            coords = pose.residue(i + 1).xyz('CA')
            data.append(" %9.3f %9.3f %9.3f" % (coords.x, coords.y, coords.z))

        if self.psize is None:
            self.psize = len(data)

        data.insert(0, "%8d %10.3f %6d %6d" % (
            self.psize, score, 1, 1)
        )

        with open(out, 'wt') as f:
            for data_line in data:
                f.write('%s\n' % data_line)

        self.tra1_filenames.append(filename)

    def init_tra1in_file(self):
        out = os.path.join(self.working_folder, 'tra.in')

        with open(out, 'wt') as f:
            f.write("%d %d %d\n" % (len(self.pop), -1, 1))

            for filename in self.tra1_filenames:
                f.write("%s\n" % filename)

    def init_rmsinp_file(self):
        out = os.path.join(self.working_folder, 'rmsinp')

        with open(out, 'wt') as f:
            f.write("%d %d\n" % (1, self.psize))
            f.write("%d\n" % (self.psize))

    def init_seq_file(self):
        out = os.path.join(self.working_folder, 'seq.dat')
        pose = self.pop[0].pose

        with open(out, 'wt') as f:
            for n in range(pose.total_residue()):
                f.write("%4d %3s\n" % (n + 1, pose.residue(n + 1).name3()))

    def init_ca_file(self):
        # TODO, leave only `ATOM` clauses, othewise Spicker crashes
        # native = self.de.rosetta_pack.native
        # out = os.path.join(self.working_folder, 'CA')
        # native.dump_pdb(out)
        pass

    def call_spicker(self):
        os.chdir(self.working_folder)
        subprocess.check_output(['../spicker'])
        os.chdir('../')

    def parse_results(self):
        n_cluster_re = re.compile('Number of clusters:')
        n_clusters = None
        cluster_summary_re = re.compile('B------------')
        filenames = []

        with open(os.path.join(self.working_folder, 'rst.dat')) as f:
            lines = f.readlines()
            for line in lines:
                if n_cluster_re.search(line):
                    n_clusters = int(line.strip().split(':')[1].strip())
                    break

            for index, line in enumerate(lines):
                if cluster_summary_re.search(line):
                    for offset in range(1, n_clusters + 1):
                        tokens = re.sub(' {2,}', ' ', lines[index + offset].strip()).split(' ')
                        filename = tokens[-1]
                        filenames.append(filename)

        for filename in filenames:
            self.centroid_index_in_pop.append(self.cluster_map[filename])

    def create_working_folder(self):
        os.mkdir(self.working_folder)

    def remove_working_folder(self):
        # rmtree(self.working_folder)
        pass


def random_string():
    char_set = string.ascii_uppercase + string.digits
    return ''.join(random.sample(char_set * 6, 6))
