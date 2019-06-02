import time
import sys


class Repacker:
    def __init__(self, de=None):
        self.de = de
        self.inject_deps()

    def inject_deps(self):
        self.rosetta_pack = self.de.rosetta_pack
        self.pop = self.de.pop

    def run_repack(self):
        mode = self.de.repack_mode

        if mode == 'all':
            self.run_repack_for_all()

        if mode == 'cluster':
            self.run_repack_for_clusters()

        if mode == 'best':
            self.run_repack_for_best()

    def run_repack_for_all(self):
        print('[REPACK] running for all')
        self.run_repack_for_list(range(self.de.pop_size))

    def run_repack_for_clusters(self):
        print('[REPACK] running for cluster centroids')
        indexes = self.de.spicker.centroid_index_in_pop
        self.run_repack_for_list(indexes)

    def run_repack_for_list(self, indexes):
        best_index = self.de.best_index

        self.dump_pbd_from_index_list(indexes)

        for index in indexes:
            print('[REPACK] running for %4d' % index)
            sys.stdout.flush()

            start_time = time.time()
            is_best = index == best_index

            individual = self.pop[index]
            rmsd = self.rosetta_pack.get_rmsd_from_native(individual.pose)
            oldscore = individual.score
            score = individual.repack()

            name_preffix = '%04d' % index
            name = self.rosetta_pack.protein_loader.original + '/' + \
                ("%s_repacked_%05d_" % (name_preffix, self.de.it)) + \
                self.de.name_suffix + ".pdb"
            individual.repacked.dump_pdb(name)

            tm_before = individual.run_tmscore()
            self.rosetta_pack.run_tmscore(name=name)
            tm_after = self.rosetta_pack.get_tmscore()

            end_time = time.time()

            data = self.build_data(tm_before, tm_after, rmsd, oldscore, score, start_time, end_time)
            self.log(
                individual,
                data,
                preffix='%04d' % index
            )

            if is_best:
                self.fake_repack_for_best(individual, data)

    def dump_pbd_from_index_list(self, indexes):
        for index in indexes:
            name = self.rosetta_pack.protein_loader.original + '/' + \
                ("best_%05d_%05d_" % (index, self.de.it)) + self.de.name_suffix + ".pdb"
            self.pop[index].pose.dump_pdb(name)

    def fake_repack_for_best(self, individual, data):
        name_preffix = 'best'
        name = self.rosetta_pack.protein_loader.original + '/' + \
            ("%s_repacked_%05d_" % (name_preffix, self.de.it)) + \
            self.de.name_suffix + ".pdb"
        individual.repacked.dump_pdb(name)

        self.log(
            individual,
            data
        )

    def run_repack_for_best(self):
        print('[REPACK] running for best')

        start_time = time.time()
        rmsd = self.rosetta_pack.get_rmsd_from_native(self.pop[self.de.best_index].pose)
        oldscore = self.de.best_score
        score = self.pop[self.de.best_index].repack()

        name = self.rosetta_pack.protein_loader.original + '/' + \
            ("best_repacked_%05d_" % self.de.it) + self.de.name_suffix + ".pdb"
        self.pop[self.de.best_index].repacked.dump_pdb(name)

        tm_before = self.pop[self.de.best_index].run_tmscore()
        self.rosetta_pack.run_tmscore(name=name)
        tm_after = self.rosetta_pack.get_tmscore()

        end_time = time.time()

        data = self.build_data(tm_before, tm_after, rmsd, oldscore, score, start_time, end_time)
        self.log(
            self.pop[self.de.best_index],
            data
        )

    def build_data(self, tm_before, tm_after, rmsd, oldscore, score, start_time, end_time):
        data = {}
        data['tm_before'] = tm_before
        data['tm_after'] = tm_after
        data['rmsd'] = rmsd
        data['oldscore'] = oldscore
        data['score'] = score
        data['start_time'] = start_time
        data['end_time'] = end_time

        return data

    def log(self, individual, data, preffix=''):
        tm_before, tm_after = data['tm_before'], data['tm_after']
        rmsd = data['rmsd']
        oldscore = data['oldscore']
        score = data['score']
        start_time, end_time = data['start_time'], data['end_time']

        base_name = "repack_" if preffix == '' else "%s_repack_" % preffix

        name = self.rosetta_pack.protein_loader.original + '/' + \
            base_name + self.de.name_suffix + ".dat"

        with open(name, 'w') as f:
            f.write('repack_time:        %12.4f\n' % (end_time - start_time))
            f.write('score:              %12.4f\n' % oldscore)
            f.write('scorefxn:           %12.4f\n' % score)
            f.write('rmsd_after:         %12.4f\n' %
                    (self.rosetta_pack.get_rmsd_from_native(individual.repacked)))
            f.write('rmsd_before:        %12.4f\n' % rmsd)
            f.write('rmsd_change:        %12.4f\n' %
                    (rmsd - self.rosetta_pack.get_rmsd_from_native(individual.repacked)))
            f.write('tm_score_before:    %12.4f\n' % tm_before['tm_score'])
            f.write('maxsub_before:      %12.4f\n' % tm_before['maxsub'])
            f.write('gdt_ts_before:      %12.4f\n' % tm_before['gdt_ts'][0])
            f.write('gdt_ha_before:      %12.4f\n' % tm_before['gdt_ha'][0])
            f.write('tm_score_after:     %12.4f\n' % tm_after['tm_score'])
            f.write('maxsub_after:       %12.4f\n' % tm_after['maxsub'])
            f.write('gdt_ts_after:       %12.4f\n' % tm_after['gdt_ts'][0])
            f.write('gdt_ha_after:       %12.4f\n' % tm_after['gdt_ha'][0])
            f.write('tm_score_change:    %12.4f\n' % (tm_before['tm_score'] - tm_after['tm_score']))
            f.write('maxsub_change:      %12.4f\n' % (tm_before['maxsub'] - tm_after['maxsub']))
            f.write('gdt_ts_change:      %12.4f\n' % (tm_before['gdt_ts'][0] - tm_after['gdt_ts'][0]))
            f.write('gdt_ha_change:      %12.4f\n' % (tm_before['gdt_ha'][0] - tm_after['gdt_ha'][0]))
            f.write(
                'gdt_ts_info_before: %12.4f %12.4f %12.4f %12.4f\n' %
                (tm_before['gdt_ts'][1][0],
                 tm_before['gdt_ts'][1][1],
                    tm_before['gdt_ts'][1][2],
                    tm_before['gdt_ts'][1][3]))
            f.write(
                'gdt_ts_info_after:  %12.4f %12.4f %12.4f %12.4f\n' %
                (tm_after['gdt_ts'][1][0],
                 tm_after['gdt_ts'][1][1],
                    tm_after['gdt_ts'][1][2],
                    tm_after['gdt_ts'][1][3]))
            f.write(
                'gdt_ha_info_before: %12.4f %12.4f %12.4f %12.4f\n' %
                (tm_before['gdt_ha'][1][0],
                 tm_before['gdt_ha'][1][1],
                    tm_before['gdt_ha'][1][2],
                    tm_before['gdt_ha'][1][3]))
            f.write(
                'gdt_ha_info_after:  %12.4f %12.4f %12.4f %12.4f\n' %
                (tm_after['gdt_ha'][1][0],
                 tm_after['gdt_ha'][1][1],
                    tm_after['gdt_ha'][1][2],
                    tm_after['gdt_ha'][1][3]))
