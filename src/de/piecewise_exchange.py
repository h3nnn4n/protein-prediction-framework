import random
import sys


# TODO: Add config
# TODO: Run experiments
# FIXME: I am pretty sure this only works with the centroid model

class PiecewiseExchange:
    def __init__(self, de=None):
        self.de = de

    def split_regions(self):
        trial = self.de.trial
        ss_pred = self.de.rosetta_pack.ss_pred
        split_points = []

        for i in range(1, len(ss_pred)):
            if ss_pred[i - 1] != ss_pred[i]:
                sidechain_size = trial.bounds.getNumSideChainAngles(ss_pred[i - 1])
                split_points.append(i * 3 + sidechain_size)

        return split_points

    def random_search(self):
        n_iters = 50

        print('starting')
        sys.stdout.flush()
        for iter_count in range(n_iters):
            print("%6d" % (iter_count), end=' ')

            # score, rmsd = self.random_piecewise_search_with_stage2(n=100)

            # score, rmsd = self.random_piecewise_search_with_stage2_then_repack(n=100)
            # scorefxn = self.de.trial.repacked.energies().total_energy()

            _, rmsd = self.random_piecewise_search_with_stage2_then_repack(n=100)
            # score, rmsd = self.random_piecewise_search_with_repack()
            scorefxn = self.de.trial.repacked.energies().total_energy()
            print("%8.3f %8.3f" % (scorefxn, rmsd))
            sys.stdout.flush()

    def scramble_pop(self):
        pop_size = self.de.pop_size

        split_regions = self.split_regions()
        split_regions = [0] + split_regions + [len(self.de.rosetta_pack.ss_pred)]

        for _ in range(pop_size):
            for i in range(len(split_regions) - 1):
                region = split_regions[i], split_regions[i + 1]

                p1, p2 = random.sample(range(pop_size), k=2)

                self.swap_region(region, p1, p2)

        self.update_poses()
        self.update_scores()

    def swap_region(self, region, p1, p2):
        start, end = region

        ind1 = self.de.pop[p1]
        ind2 = self.de.pop[p2]

        t = ind1.angles[start:end].copy()
        ind1.angles[start:end] = ind2.angles[start:end]
        ind2.angles[start:end] = t

    def random_piecewise_search_with_repack(self):
        rp = self.de.rosetta_pack
        self.random_piecewise_search()
        score = self.de.trial.repack()
        rmsd = rp.get_rmsd_from_native(self.de.trial.pose)

        return score, rmsd

    def random_piecewise_search_with_stage2(self, n=25, temp=1.5, mode='3s'):
        rp = self.de.rosetta_pack
        self.random_piecewise_search()

        self.stage2()

        score = self.de.trial.score
        rmsd = rp.get_rmsd_from_native(self.de.trial.pose)

        return score, rmsd

    def random_piecewise_search_with_stage2_then_repack(self, n=25, temp=1.5, mode='3s'):
        rp = self.de.rosetta_pack
        self.random_piecewise_search()

        self.stage2()

        self.de.trial.repack()
        self.de.trial.eval()
        score = self.de.trial.score
        rmsd = rp.get_rmsd_from_native(self.de.trial.pose)

        return score, rmsd

    def random_piecewise_search(self):
        trial = self.de.trial
        pop_size = self.de.pop_size

        split_regions = self.split_regions()
        split_regions = [0] + split_regions + [len(self.de.rosetta_pack.ss_pred)]

        for i in range(len(split_regions) - 1):
            region = split_regions[i], split_regions[i + 1]
            p1, _ = random.sample(range(pop_size), k=2)
            self.add_region(region, p1)

        trial.update_pose_from_angles()
        trial.eval()

    def add_region(self, region, p1):
        start, end = region

        ind1 = self.de.pop[p1]
        trial = self.de.trial

        trial.angles[start:end] = ind1.angles[start:end]

    def update_poses(self):
        for p in self.de.pop:
            p.update_pose_from_angles()

    def update_scores(self):
        for p in self.de.pop:
            p.eval()

    def stage2(self):
        trial = self.de.trial
        pd = self.de.rosetta_pack

        score = pd.get_score_function(trial.score_function_name)
        temp = 1.5
        n = 250

        mc = pd.get_new_mc(trial.pose, score, temp)
        mover = pd.get_3mer_smooth()
        trial_mover = pd.get_new_trial_mover(mover, mc)
        rep_mover = pd.get_new_rep_mover(trial_mover, n)

        mc.set_temperature(temp)

        trial.eval()
        score_before = trial.score
        rep_mover.apply(trial.pose)

        trial.pose.assign(mc.lowest_score_pose())
        trial.update_angle_from_pose()
        trial.eval()

        score_after = trial.score

        print("%8.3f %8.3f" % (score_before, score_after), end=' ')
