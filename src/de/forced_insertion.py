from random import uniform

# TODO: Test coverage
# TODO: Check for all atom


class ForcedInsertion:
    def __init__(self, de=None, mode=None, chance=0.0):
        # print("forced_frag init mode: %s  chance: %f" % (mode, chance))
        self.de = de
        self.mode = mode
        self.chance = chance

        self.frag_inserter = None

        self.initialize_rosetta()

    def run(self):
        if not self.can_run():
            return

        de = self.de

        for p in de.pop:
            chance = uniform(0, 1)

            if chance > self.chance:
                continue

            self.apply_insertion(p)

    def can_run(self):
        return self.de is not None and self.mode is not None and self.chance > 0

    def apply_insertion(self, individual):
        # TODO add some metrics and logging here

        # score_before = individual.score

        self.frag_inserter.apply(individual.pose)
        individual.update_angle_from_pose()
        individual.eval()

        # score_after = individual.score

        # print("forced frag: %8.3f %8.3f %8.3f" % (
        #     score_before,
        #     score_after,
        #     score_after - score_before
        # ))

    def initialize_rosetta(self):
        if not self.can_run():
            return

        self.frag_inserter = self.get_frag_mover()

    def get_frag_mover(self):
        rp = self.de.rosetta_pack

        if self.mode == '9mer' or self.mode == 'frag9':
            return rp.get_9mer()

        raise ValueError(
            "%s is not a valid frag mover or it is not supported it" % self.mode
        )
