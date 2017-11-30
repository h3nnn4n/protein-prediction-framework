import numpy as np


class ProteinData:
    def __init__(self, rosetta_pack, allatom=False):
        self.rosetta_pack = rosetta_pack
        self.pose = rosetta_pack.get_new_pose()
        self.bounds = rosetta_pack.bounds
        self.allatom = allatom

        if allatom:
            self.rosetta_pack.get_allatom_switch().apply(self.pose)

        self.nsca = sum([(self.bounds.getNumSideChainAngles(ss) + 3) for ss in self.rosetta_pack.target])
        self.angles = np.zeros(self.nsca)
        self.score = None

        index = 0
        for k, ss in enumerate(rosetta_pack.ss_pred):
            n_sidechain_angles = self.bounds.getNumSideChainAngles(rosetta_pack.target[k])
            phi, psi, omega = self.bounds.generateRandomAngles(ss)

            # print(index, k, n_sidechain_angles, self.nsca)

            self.pose.set_phi(k + 1, phi)
            self.pose.set_psi(k + 1, psi)
            self.pose.set_omega(k + 1, omega)

            self.angles[index + 0] = phi
            self.angles[index + 1] = psi
            self.angles[index + 2] = omega

            aa = self.rosetta_pack.target[k]
            angles = self.bounds.generateRandomSidechainAngles(aa)
            for kk, vv in enumerate(angles):
                if self.allatom and n_sidechain_angles > 0:
                    self.angles[index + 3 + kk] = vv
                    self.pose.set_chi(kk + 1, k + 1, vv)

            print("%4d %8.3f %8.3f %8.3f " % (k, self.angles[index + 0], self.angles[index + 1], self.angles[index + 2]), end='')
            for i in range(len(angles)):
                print(" %8.3f" % self.angles[index + 3 + i], end='')
            print()

            index += 3 + n_sidechain_angles

        self.fix_bounds()
        self.eval()
        print()

        # print('Finished INIT', self.allatom)

    def reset(self):
        index = 0
        for k, ss in enumerate(self.rosetta_pack.ss_pred):
            n_sidechain_angles = self.bounds.getNumSideChainAngles(self.rosetta_pack.target[k])
            phi, psi, omega = self.bounds.generateRandomAngles(ss)

            self.pose.set_phi(k + 1, phi)
            self.pose.set_psi(k + 1, psi)
            self.pose.set_omega(k + 1, omega)

            self.angles[index + 0] = phi
            self.angles[index + 1] = psi
            self.angles[index + 2] = omega

            if self.allatom and n_sidechain_angles > 0:
                aa = self.rosetta_pack.target[k]
                angles = self.bounds.generateRandomSidechainAngles(aa)
                for kk, vv in enumerate(angles):
                    self.angles[index + 3 + kk] = vv
                    self.pose.set_chi(kk + 1, k + 1, vv)

            index += 3 + n_sidechain_angles

        self.eval()

        self.fix_bounds()

    def new_angles(self, angles):
        index = 0
        for k, ss in enumerate(self.rosetta_pack.ss_pred):
            n_sidechain_angles = self.bounds.getNumSideChainAngles(self.rosetta_pack.target[k])
            phi, psi, omega = angles[index + 0], angles[index + 1], angles[index + 2]

            self.pose.set_phi(k + 1, phi)
            self.pose.set_psi(k + 1, psi)
            self.pose.set_omega(k + 1, omega)

            self.angles[index + 0] = phi
            self.angles[index + 1] = psi
            self.angles[index + 2] = omega

            if self.allatom and n_sidechain_angles > 0:
                for j in range(n_sidechain_angles):
                    vv = angles[index + 3 + j]
                    self.angles[index + 3 + j] = vv
                    self.pose.set_chi(j + 1, k + 1, vv)

            index += 3 + n_sidechain_angles

    def eval(self, score=None):
        if score is None and self.pose.is_centroid():
            # print('centroid')
            score = self.rosetta_pack.get_score3()

        if score is None and not self.pose.is_centroid():
            # print('allatom')
            score = self.rosetta_pack.get_scorefxn()

        self.score = score(self.pose)
        # self.rosetta_pack.pymover.apply(self.pose)

    def fix_bounds(self):
        index = 0
        for n, ss in enumerate(self.rosetta_pack.ss_pred):
            phi = self.angles[index + 0]
            psi = self.angles[index + 1]
            # omega = self.angles[index + 2]
            n_sidechain_angles = self.bounds.getNumSideChainAngles(self.rosetta_pack.target[n])

            if phi < self.bounds.getSecondaryLowerBound(ss, 0):
                phi = self.bounds.getSecondaryLowerBound(ss, 0)
                self.pose.set_phi(n + 1, phi)
                self.angles[index + 0] = phi

            if phi > self.bounds.getSecondaryUpperBound(ss, 0):
                phi = self.bounds.getSecondaryUpperBound(ss, 0)
                self.pose.set_phi(n + 1, phi)
                self.angles[index + 0] = phi

            if psi < self.bounds.getSecondaryLowerBound(ss, 1):
                psi = self.bounds.getSecondaryLowerBound(ss, 1)
                self.pose.set_psi(n + 1, psi)
                self.angles[index + 1] = psi

            if psi > self.bounds.getSecondaryUpperBound(ss, 1):
                psi = self.bounds.getSecondaryUpperBound(ss, 1)
                self.pose.set_psi(n + 1, psi)
                self.angles[index + 1] = psi

            index += 3 + n_sidechain_angles

    def update_angle_from_pose(self):
        # for n, _ in enumerate(self.rosetta_pack.ss_pred):
            # phi, psi, omega = self.pose.phi(n + 1), self.pose.psi(n + 1), self.pose.omega(n + 1)

#             if self.angles[n * 3 + 0] - phi > 0.001 or self.angles[n * 3 + 1] - psi > 0.001:
#                 print("%d" % n)
#                 print("%8.3f %8.3f %8.3f" % (phi, psi, omega))
#                 print("%8.3f %8.3f %8.3f" % (self.angles[n * 3 + 0], self.angles[n * 3 + 1], self.angles[n * 3 + 2]))

            # self.angles[n * 3 + 0] = phi
            # self.angles[n * 3 + 1] = psi
            # self.angles[n * 3 + 2] = omega

        index = 0
        for k, ss in enumerate(self.rosetta_pack.ss_pred):
            n_sidechain_angles = self.bounds.getNumSideChainAngles(self.rosetta_pack.target[k])
            phi, psi, omega = self.pose.phi(k + 1), self.pose.psi(k + 1), self.pose.omega(k + 1)

            self.angles[index + 0] = phi
            self.angles[index + 1] = psi
            self.angles[index + 2] = omega

            if self.allatom and n_sidechain_angles > 0:
                for j in range(n_sidechain_angles):
                    vv = self.pose.chi(j + 1, k + 1)
                    self.angles[index + 3 + j] = vv

            index += 3 + n_sidechain_angles

    def stage1_mc(self, n=100, temp=2.0):
        allatom = self.allatom

        pd = self.rosetta_pack
        pd.set_starting_pose(self.angles)

        best = self.pose
        score = pd.get_score0()

        r = None
        if allatom:
            r = pd.get_sidechain_recover()(best)
            pd.get_centroid_switch().apply(best)

        mc = pd.get_new_mc(best, score, temp)

        seq = pd.get_new_seq_mover()
        seq.add_mover(pd.get_9mer())
        # seq.add_mover(pd.get_9mer_smooth())
        # seq.add_mover(pd.get_3mer())
        # seq.add_mover(pd.get_3mer_smooth())

        trial = pd.get_new_trial_mover(seq, mc)

        folding = pd.get_new_rep_mover(trial, n)
        folding.apply(best)

        mc.recover_low(best)
        mc.reset(best)

        if allatom:
            pd.get_allatom_switch().apply(best)
            r.apply(best)
            pd.get_packer().apply(best)

    def stage2_mc(self, n=100, temp=2.0):
        allatom = self.allatom

        pd = self.rosetta_pack
        pd.set_starting_pose(self.angles)

        best = self.pose
        score = pd.get_score3()

        r = None
        if allatom:
            r = pd.get_sidechain_recover()(best)
            pd.get_centroid_switch().apply(best)

        mc = pd.get_new_mc(best, score, temp)

        seq = pd.get_new_seq_mover()
        # seq.add_mover(pd.get_9mer())
        # seq.add_mover(pd.get_9mer_smooth())
        # seq.add_mover(pd.get_3mer())
        seq.add_mover(pd.get_3mer_smooth())

        trial = pd.get_new_trial_mover(seq, mc)

        folding = pd.get_new_rep_mover(trial, n)
        folding.apply(best)

        mc.recover_low(best)
        mc.reset(best)

        if allatom:
            pd.get_allatom_switch().apply(best)
            r.apply(best)
            pd.get_packer().apply(best)
