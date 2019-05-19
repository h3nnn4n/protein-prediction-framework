import time


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

        if mode == 'best':
            self.run_repack_for_best()

    def run_repack_for_all(self):
        pass

    def run_repack_for_best(self):
        rmsd = self.rosetta_pack.get_rmsd_from_native(self.pop[self.de.best_index].pose)
        oldscore = self.de.best_score
        score = self.pop[self.de.best_index].repack()

        name = self.rosetta_pack.protein_loader.original + '/' + \
            ("best_repacked_%05d_" % self.de.it) + self.de.name_suffix + ".pdb"
        self.pop[self.de.best_index].repacked.dump_pdb(name)

        repack_name = name

        tm_before = self.pop[self.de.best_index].run_tmscore()
        self.rosetta_pack.run_tmscore(name=repack_name)
        tm_after = self.rosetta_pack.get_tmscore()

        self.repack_time = time.time()
        name = self.rosetta_pack.protein_loader.original + '/' + "repack_" + self.de.name_suffix + ".dat"
        with open(name, 'w') as f:
            f.write('repack_time:        %12.4f\n' % (self.repack_time - self.de.end_time))
            f.write('score:              %12.4f\n' % oldscore)
            f.write('scorefxn:           %12.4f\n' % score)
            f.write('rmsd_after:         %12.4f\n' %
                    (self.rosetta_pack.get_rmsd_from_native(self.pop[self.de.best_index].repacked)))
            f.write('rmsd_before:        %12.4f\n' % rmsd)
            f.write('rmsd_change:        %12.4f\n' %
                    (rmsd - self.rosetta_pack.get_rmsd_from_native(self.pop[self.de.best_index].repacked)))
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
