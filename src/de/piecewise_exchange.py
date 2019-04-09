import random


# TODO: Add config
# TODO: Random search mode
# TODO: Random Search with
# TODO: Run experiments

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

        # TODO: Use slices
        for i in range(start, end):
            t = ind1.angles[i]
            ind1.angles[i] = ind2.angles[i]
            ind2.angles[i] = t

    def update_poses(self):
        for p in self.de.pop:
            p.update_pose_from_angles()

    def update_scores(self):
        for p in self.de.pop:
            p.eval()
