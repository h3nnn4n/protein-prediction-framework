import sys
import time
import random

from shutil import copyfile

from local_search.hooke_jeeves import hooke


class HookeJeevesPostprocessing:
    def __init__(self, de=None):
        self.de = de
        self.inject_deps()

    def inject_deps(self):
        self.pop = self.de.pop
        self.rosetta_pack = self.de.rosetta_pack

    def run_hooke_jeeves(self):
        if not self.de.hooke_jeeves_postprocessing:
            return

        mode = self.de.hooke_jeeves_postprocessing_mode

        if mode == 'all':
            self.run_hooke_jeeves_for_all()

        if mode == 'best':
            self.run_hooke_jeeves_for_best()

    def run_hooke_jeeves_for_all(self):
        for index, individual in enumerate(self.de.pop):
            print('running hooke_jeeves for %4d' % index)
            sys.stdout.flush()

            self.run_hooke_jeeves_for_individual(
                individual,
                preffix=('%04d' % index),
                is_best=(index == self.de.best_index)
            )

    def run_hooke_jeeves_for_best(self):
        self.run_hooke_jeeves_for_individual(self.pop[self.de.best_index], preffix='best')

    def run_hooke_jeeves_for_individual(self, individual, preffix='', is_best=False):
        name = self.rosetta_pack.protein_loader.original + '/' + \
            ("hooke-jeeves_%s_" % preffix) + self.de.name_suffix + ".dat"
        hooke_start_time = time.time()

        score_before, score_after = individual.score, None
        rmsd_before, rmsd_after = self.rosetta_pack.get_rmsd_from_native(individual.pose), None

        scores = [score_before]
        rmsds = [rmsd_before]
        spent_evals = [0]

        while True:
            old_score = individual.score
            evals, _ = hooke(
                individual,
                eps=1e-04,
                rho=random.random(),
                itermax=250
            )
            new_score = individual.score
            rmsd = self.rosetta_pack.get_rmsd_from_native(individual.pose)

            scores.append(new_score)
            spent_evals.append(evals)
            rmsds.append(rmsd)

            if old_score - new_score < 0.01:
                break

        score_after = individual.score
        rmsd_after = self.rosetta_pack.get_rmsd_from_native(individual.pose)

        hooke_end_time = time.time()

        with open(name, 'w') as f:
            f.write('hooke_time:        %12.4f\n' % (hooke_end_time - hooke_start_time))
            f.write('score_before:      %12.4f\n' % score_before)
            f.write('score_after:       %12.4f\n' % score_after)
            f.write('rmsd_before:       %12.4f\n' % rmsd_before)
            f.write('rmsd_after:        %12.4f\n' % rmsd_after)
            f.write('spent_evals:       %12d\n' % sum(spent_evals))

            for score_index, score in enumerate(scores):
                f.write('score_%02d:          %12.4f\n' % (score_index + 1, score))

            for rmsd_index, rmsd in enumerate(rmsds):
                f.write('rmsd_%02d:           %12.4f\n' % (rmsd_index + 1, rmsd))

            for eval_index, evals in enumerate(spent_evals):
                f.write('spent_evals_%02d:    %12d\n' % (eval_index + 1, evals))

        if is_best:
            src = name
            dest = self.rosetta_pack.protein_loader.original + '/' + \
                "hooke-jeeves_best_" + self.de.name_suffix + ".dat"
            copyfile(src, dest)
