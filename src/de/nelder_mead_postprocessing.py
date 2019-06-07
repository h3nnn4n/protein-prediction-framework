import sys
import time
import random

from shutil import copyfile

from local_search.nelder_mead import nelder_mead


class NelderMeadPostprocessing:
    def __init__(self, de=None):
        self.de = de
        self.inject_deps()

    def inject_deps(self):
        self.pop = self.de.pop
        self.rosetta_pack = self.de.rosetta_pack

    def run_nelder_mead(self):
        if not self.de.nelder_mead_postprocessing:
            return

        mode = self.de.nelder_mead_postprocessing_mode

        if mode == 'all':
            self.run_nelder_mead_for_all()

        if mode == 'cluster':
            self.run_nelder_mead_for_clusters()

        if mode == 'best':
            self.run_nelder_mead_for_best()

    def run_nelder_mead_for_all(self):
        print('[NELDER MEAD] running for all')

        for index, individual in enumerate(self.de.pop):
            print('[NELDER MEAD] running for %4d' % index)
            sys.stdout.flush()

            self.run_nelder_mead_for_individual(
                individual,
                preffix=('%04d' % index),
                is_best=(index == self.de.best_index)
            )

    def run_nelder_mead_for_clusters(self):
        print('[NELDER MEAD] running for cluster centroids')

        indexes = self.de.spicker.centroid_index_in_pop
        for index in indexes:
            print('[NELDER MEAD] running for %4d' % index)
            sys.stdout.flush()
            individual = self.pop[index]

            self.run_nelder_mead_for_individual(
                individual,
                preffix=('%04d' % index),
                is_best=(index == self.de.best_index)
            )

    def run_nelder_mead_for_best(self):
        print('[NELDER MEAD] running for best')

        self.run_nelder_mead_for_individual(self.pop[self.de.best_index], preffix='best')

    def run_nelder_mead_for_individual(self, individual, preffix='', is_best=False):
        name = self.rosetta_pack.protein_loader.original + '/' + \
            ("nelder-mead_%s_" % preffix) + self.de.name_suffix + ".dat"
        nelder_start_time = time.time()

        score_before, score_after = individual.score, None
        rmsd_before, rmsd_after = self.rosetta_pack.get_rmsd_from_native(individual.pose), None

        scores = [score_before]
        rmsds = [rmsd_before]
        spent_evals = [0]

        while True:
            old_score = individual.score
            evals, _ = nelder_mead(
                individual,
                eps=1e-04,
                max_evals=500
            )
            new_score = individual.score
            rmsd = self.rosetta_pack.get_rmsd_from_native(individual.pose)

            scores.append(new_score)
            spent_evals.append(evals)
            rmsds.append(rmsd)

            if old_score - new_score < 0.1:
                break

        score_after = individual.score
        rmsd_after = self.rosetta_pack.get_rmsd_from_native(individual.pose)

        nelder_end_time = time.time()

        with open(name, 'w') as f:
            f.write('nelder_time:        %12.4f\n' % (nelder_end_time - nelder_start_time))
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
                "nelder-mead_best_" + self.de.name_suffix + ".dat"
            copyfile(src, dest)
