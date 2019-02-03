import time
import sys
import datetime
import string
import pyrosetta
import random
from rosetta_pack import RosettaPack


class ClassicAbinitio:
    def __init__(self, pname='1zdd'):
        self.pname = pname
        self.cname = '%s_classic-abinitio' % self.pname

        self.rp = RosettaPack(self.pname)
        self.score3 = self.rp.get_score3()
        self.pose = self.rp.get_new_pose()
        self.abinitio = self.rp.get_rosetta_abinitio_protocol()
        self.repacked = None
        self.rosetta_pack = self.rp

        self.setup_logs()

    def setup_logs(self):
        self.init_time = time.time()
        self.now = datetime.datetime.now()
        self.start_time = self.init_time

        char_set = string.ascii_uppercase + string.digits
        r_string = ''.join(random.sample(char_set * 6, 6))

        now = self.now
        self.name_suffix = "_%s__%s__%04d_%02d_%02d__%02d_%02d_%02d__%s" % (
            self.pname,
            self.cname,
            now.year,
            now.month,
            now.day,
            now.hour,
            now.minute,
            now.second,
            r_string
        )

        self.stats_name = self.rosetta_pack.protein_loader.original + '/' + "stats_" + self.name_suffix + ".dat"
        with open(self.stats_name, 'wt') as f:
            f.flush()

        self.config_name = self.rosetta_pack.protein_loader.original + '/' + "parameters_" + self.name_suffix + ".yaml"
        with open(self.config_name, 'w') as f:
            f.write('hello world\n')
            f.flush()

    def run(self, factor=1):
        self.log()
        self.abinitio.set_cycles(factor)
        self.abinitio.apply(self.pose)
        self.log(1)
        self.end_time = time.time()
        self.dump_pbd_best()
        self.dump_pdb_repacked()

    def dump_pbd_best(self, it=1):
        name = self.rosetta_pack.protein_loader.original + '/' + ("best_%05d_" % it) + self.name_suffix + ".pdb"
        self.pose.dump_pdb(name)

    def dump_pdb_repacked(self, it=1):
        rmsd = self.rosetta_pack.get_rmsd_from_native(self.pose)
        oldscore = self.score3(self.pose)

        score = self.repack()

        name = self.rosetta_pack.protein_loader.original + '/' + ("best_repacked_%05d_" % it) + self.name_suffix + ".pdb"
        self.repacked.dump_pdb(name)

        repack_name = name

        self.rosetta_pack.pose.assign(self.pose)
        self.rosetta_pack.run_tmscore(pose=self.repacked)
        tm_before = self.rosetta_pack.get_tmscore()

        self.rosetta_pack.run_tmscore(name=repack_name)
        tm_after = self.rosetta_pack.get_tmscore()

        self.repack_time = time.time()
        name = self.rosetta_pack.protein_loader.original + '/' + "repack_" + self.name_suffix + ".dat"
        with open(name, 'w') as f:
            f.write('repack_time:        %12.4f\n' % (self.repack_time - self.end_time))
            f.write('score:              %12.4f\n' % oldscore)
            f.write('scorefxn:           %12.4f\n' % score)
            f.write('rmsd_after:         %12.4f\n' % (self.rosetta_pack.get_rmsd_from_native(self.repacked)))
            f.write('rmsd_before:        %12.4f\n' % rmsd)
            f.write('rmsd_change:        %12.4f\n' % (rmsd - self.rosetta_pack.get_rmsd_from_native(self.repacked)))
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
            f.write('gdt_ts_info_before: %12.4f %12.4f %12.4f %12.4f\n' % (tm_before['gdt_ts'][1][0], tm_before['gdt_ts'][1][0], tm_before['gdt_ts'][1][0], tm_before['gdt_ts'][1][0]))
            f.write('gdt_ha_info_after:  %12.4f %12.4f %12.4f %12.4f\n' % (tm_after['gdt_ha'][1][0], tm_after['gdt_ha'][1][0], tm_after['gdt_ha'][1][0], tm_after['gdt_ha'][1][0]))
            f.flush()

    def log(self, it=0):
        rmsd = self.rp.get_rmsd_from_native(self.pose)
        energy = self.score3(self.pose)
        evals = self.abinitio.total_trials()

        if it == 0:
            evals = 0

        data = [
            ('%8d', evals),
            ('%8d', it),
            ('%8.4f', energy),
            ('%8.4f', energy),
            ('%8.4f', 0),
            ('%8.4f', rmsd),
            ('%8.4f', rmsd),
            ('%8.4f', 0),
            ('%10.2f', (time.time() - self.start_time)),
            ('%8.5f', 0),
            ('%8.2f', 0),
            ('%8.5f', 0),
            ('%8.2f', 0),
        ]

        string = ''
        for k, (a, b) in enumerate(data):
            try:
                string += a % b
                string += ' '
            except Exception:
                print(k, a, b)

        print(string)

        with open(self.stats_name, 'at') as f:
            f.write(string + '\n')
            f.flush()

    def repack(self):
        repack = self.rp.get_fast_relax()
        best = self.rp.copy_pose_to_allatom(self.pose)
        repack.apply(best)
        self.repacked = best
        return self.rp.get_scorefxn()(best)


if __name__ == '__main__':
    pname = '1zdd'
    factor = 1
    if len(sys.argv) > 1:
        pname = sys.argv[1]

    if len(sys.argv) > 2:
        factor = int(sys.argv[2])

    c = ClassicAbinitio(pname)
    c.run(factor)
