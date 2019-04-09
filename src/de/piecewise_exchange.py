import random


# TODO: Add config
# TODO: Run experiments
# FIXME: I am pretty sure this only works with the centroid model

class PiecewiseExchange:
    def __init__(self, de=None):
        self.de = de

    def split_regions(self):
        ss_pred = self.de.rosetta_pack.ss_pred
        split_points = []

        for i in range(1, len(ss_pred)):
            if ss_pred[i - 1] != ss_pred[i]:
                split_points.append(i)

        return split_points

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

        self.de.trial.stage2_mc(n=n, temp=temp, mode=mode)

        score = self.de.trial.score
        rmsd = rp.get_rmsd_from_native(self.de.trial.pose)

        return score, rmsd

    def random_piecewise_search_with_stage2_then_repack(self, n=25, temp=1.5, mode='3s'):
        rp = self.de.rosetta_pack
        self.random_piecewise_search()

        self.de.trial.stage2_mc(n=n, temp=temp, mode=mode)

        self.de.trial.repack()
        self.de.trial.eval()
        score = self.de.trial.score
        rmsd = rp.get_rmsd_from_native(self.de.trial.pose)

        return score, rmsd

    def random_piecewise_search(self):
        pop_size = self.de.pop_size

        split_regions = self.split_regions()
        split_regions = [0] + split_regions + [len(self.de.rosetta_pack.ss_pred)]

        for i in range(len(split_regions) - 1):
            region = split_regions[i], split_regions[i + 1]
            p1, _ = random.sample(range(pop_size), k=2)
            self.add_region(region, p1)

        self.de.trial.update_pose_from_angles()
        self.de.trial.eval()

    def add_region(self, region, p1):
        start, end = region

        ind1 = self.de.pop[p1]
        trial = self.de.trial

        print(start, end)

        trial.angles[start:end] = ind1.angles[start:end]

    def update_poses(self):
        for p in self.de.pop:
            p.update_pose_from_angles()

    def update_scores(self):
        for p in self.de.pop:
            p.eval()
