import numpy as np


class ProteinData:
    def __init__(self, rosetta_pack, allatom=False):
        self.rosetta_pack = rosetta_pack
        self.pose = rosetta_pack.get_new_pose()
        self.bounds = rosetta_pack.bounds
        self.allatom = allatom

        self.enable_remc = True

        self.repacked = None

        self.score_function = None

        self.mover3 = None
        self.mover3s = None
        self.mover9 = None
        self.mover9s = None
        self.moversmall = None
        self.movershear = None

        self.mc3 = None
        self.mc3s = None
        self.mc9 = None
        self.mc9s = None
        self.mcsmall = None
        self.mcshear = None

        # self.seq3 = None
        # self.seq3s = None
        # self.seq9 = None
        # self.seq9s = None
        # self.seqsmall = None
        # self.seqshear = None

        self.trial3 = None
        self.trial3s = None
        self.trial9 = None
        self.trial9s = None
        self.trialsmall = None
        self.trialshear = None

        if allatom:
            self.rosetta_pack.get_allatom_switch().apply(self.pose)

        self.nsca = sum([(self.bounds.getNumSideChainAngles(ss) + 3) for ss in self.rosetta_pack.target])
        self.angles = np.zeros(self.nsca)
        self.score = None

        self.set_score_function()
        self.reset()

        # self.fix_bounds()
        # self.eval()
        # self.print_angles()
        # print('Finished INIT', self.allatom)

    def __call__(self, angles):
        self.new_angles(angles)
        self.eval()
        return self.score

    def set_score_function(self, score='score3'):
        self.score_function_name = score
        self.score_function = self.rosetta_pack.get_score_function(score)

    def print_angles(self):
        index = 0
        for k, ss in enumerate(self.rosetta_pack.ss_pred):
            n_sidechain_angles = self.bounds.getNumSideChainAngles(self.rosetta_pack.target[k])
            print("%4d %4d %8.3f %8.3f %8.3f " % (k, index, self.angles[index + 0], self.angles[index + 1], self.angles[index + 2]), end='')
            for i in range(n_sidechain_angles):
                print(" %8.3f" % self.angles[index + 3 + i], end='')
            index += 3 + n_sidechain_angles
            print()
        print()

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

    def new_angles(self, angles):
        index = 0
        for k, _ in enumerate(self.rosetta_pack.ss_pred):
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
        self.score = self.score_function(self.pose)
        # if score is None and self.pose.is_centroid():
        #     print('centroid')
        #     score = self.rosetta_pack.get_score3()

        # if score is None and not self.pose.is_centroid():
        #     print('allatom')
        #     score = self.rosetta_pack.get_scorefxn()

        # self.score = score(self.pose)
        # self.rosetta_pack.pymover.apply(self.pose)

    def fix_bounds(self):
        return

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
        index = 0
        for k, _ in enumerate(self.rosetta_pack.ss_pred):
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
        pd = self.rosetta_pack

        score = pd.get_score0()

        mc = pd.get_new_mc(self.pose, score, temp)

        seq = pd.get_new_seq_mover()
        seq.add_mover(pd.get_9mer())
        trial = pd.get_new_trial_mover(seq, mc)
        folding = pd.get_new_rep_mover(trial, n)

        folding.apply(self.pose)

        self.pose.assign(mc.lowest_score_pose())
        self.update_angle_from_pose()
        self.eval()

    def stage2_mc(self, n=100, temp=2.0, mode='3s'):
        allatom = self.allatom

        pd = self.rosetta_pack

        score = pd.get_score_function(self.score_function_name)

        r = None
        if allatom:
            r = pd.get_sidechain_recover()(self.pose)
            pd.get_centroid_switch().apply(self.pose)

        if mode == '3':
            if self.mc3 is None:
                self.mc3 = pd.get_new_mc(self.pose, score, temp)
                self.mover3 = pd.get_3mer()
                self.trial3 = pd.get_new_trial_mover(self.mover3, self.mc3)
        elif mode == '3s':
            if self.mc3s is None:
                self.mc3s = pd.get_new_mc(self.pose, score, temp)
                self.mover3s = pd.get_3mer_smooth()
                self.trial3s = pd.get_new_trial_mover(self.mover3s, self.mc3s)
        elif mode == '9':
            if self.mc9 is None:
                self.mc9 = pd.get_new_mc(self.pose, score, temp)
                self.mover9 = pd.get_9mer()
                self.trial9 = pd.get_new_trial_mover(self.mover9, self.mc9)
        elif mode == '9s':
            if self.mc9s is None:
                self.mc9s = pd.get_new_mc(self.pose, score, temp)
                self.mover9s = pd.get_9mer_smooth()
                self.trial9s = pd.get_new_trial_mover(self.mover9s, self.mc9s)
        elif mode == 'small':
            if self.mcsmall is None:
                self.mcsmall = pd.get_new_mc(self.pose, score, temp)
                self.moversmall = pd.get_new_small_mover()
                self.trialsmall = pd.get_new_trial_mover(self.moversmall, self.mcsmall)
        elif mode == 'shear':
            if self.mcshear is None:
                self.mcshear = pd.get_new_mc(self.pose, score, temp)
                self.movershear = pd.get_new_shear_mover()
                self.trialshear = pd.get_new_trial_mover(self.movershear, self.mcshear)

        mover = None
        mc = None
        evals = 0

        if mode == '3':
            mc = self.mc3
            # mover = self.seq3
            mover = self.trial3
        elif mode == '3s':
            mc = self.mc3s
            # mover = self.seq3s
            mover = self.trial3s
        elif mode == '9':
            mc = self.mc9
            # mover = self.seq9
            mover = self.trial9
        elif mode == '9s':
            mc = self.mc9s
            # mover = self.seq9s
            mover = self.trial9s
        elif mode == 'small':
            mc = self.mcsmall
            # mover = self.seqsmall
            mover = self.trialsmall
        elif mode == 'shear':
            mc = self.mcshear
            # mover = self.seqshear
            mover = self.trialshear

        mc.set_temperature(temp)
        if not self.enable_remc:
            mc.reset(self.pose)

        self.eval()
        original = self.score
        one_more = original is None

        for _ in range(n):
            mover.apply(self.pose)
            evals += 1
            if self.pose.energies().total_energy() < original:
                if not one_more:
                    break
                else:
                    one_more = False

        self.pose.assign(mc.lowest_score_pose())
        self.score = mc.lowest_score()

        if allatom:
            pd.get_allatom_switch().apply(self.pose)
            r.apply(self.pose)
            pd.get_packer().apply(self.pose)

        self.update_angle_from_pose()

        return evals

    def repack(self):
        pd = self.rosetta_pack
        repack = pd.get_fast_relax()
        best = pd.convert_to_allatom_pose(self.pose)
        repack.apply(best)

        self.repacked = best

        return pd.get_scorefxn()(best)

    def run_tmscore(self):
        pd = self.rosetta_pack
        pd.run_tmscore()
        self.tmscore = pd.get_tmscore()
        return self.tmscore

    def copy(self, original):
        self.new_angles(original.angles)
        self.score = original.score

