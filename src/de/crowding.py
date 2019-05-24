import random


class Crowding:
    def __init__(self, de=None):
        self.de = de
        self.run_crowding = self.de.run_crowding
        self.crowding_factor = self.de.crowding_factor
        self.crowding_mode = self.de.crowding_mode

        self.safety_checks = False

        self.mean = None
        self.best_score = float('inf')

        self.inject_deps()

    def inject_deps(self):
        self.pop = self.de.pop
        self.pop_size = self.de.pop_size

    def selection(self, _):
        if self.safety_checks:
            self.update_mean()

        candidates = self.get_random_individuals()
        candidate = self.find_nearest_candidate_to_trial(candidates)

        if self.de.trial.score < candidate.score:
            self.update_sade_success()
            self.swap_candidate_and_trial(candidate)
        else:
            self.update_sade_fail()

        if self.safety_checks:
            self.check_mean()

    def swap_candidate_and_trial(self, candidate):
        swapped = False

        for k in range(self.pop_size):
            if candidate is self.pop[k]:
                t = self.pop[k]
                self.pop[k] = self.de.trial
                self.de.trial = t

                if not swapped:
                    swapped = True
                else:
                    raise Exception('Found repeated elements while doing selection!')

        assert swapped

    def find_nearest_candidate_to_trial(self, candidates):
        return self.sort_candidates_by_trial_rmsd(candidates)[0]

    def sort_candidates_by_trial_rmsd(self, candidates):
        rp = self.de.rosetta_pack
        trial_pose = self.de.trial.pose

        def sorting_key(ind):
            return rp.get_rmsd_from_pose(ind.pose, trial_pose)

        return sorted(candidates, key=sorting_key)

    def get_random_individuals(self, n=None):
        if n is None:
            n = self.crowding_factor

        return random.sample(self.de.pop, n)

    def update_sade_fail(self):
        if not self.de.sade_run:
            return

        sade_k = self.de.sade_k

        _, _ = self.de.get_f_cr()
        ind = self.de.it % self.de.sade_lp

        if self.de.sade_run:
            self.de.sade_failure_memory[ind][sade_k] += 1

    def update_sade_success(self):
        if not self.de.sade_run:
            return

        sade_k = self.de.sade_k

        _, cr = self.de.get_f_cr()
        ind = self.de.it % self.de.sade_lp

        if self.de.sade_run:
            self.de.sade_cr_memory[sade_k].append(cr)
            self.de.sade_success_memory[ind][sade_k] += 1

    def check_mean(self):
        old_mean = self.mean
        old_best = self.best_score

        self.update_mean()

        if old_mean is None:
            return

        assert self.mean <= old_mean
        assert self.best_score <= old_best

    def update_mean(self):
        self.mean = 0
        self.best_score = float('inf')

        for i in range(self.pop_size):
            self.mean += self.pop[i].score / self.pop_size
            if self.best_score is None or self.pop[i].score < self.best_score:
                self.best_score = self.pop[i].score
                self.best_index = i

        return self.mean
