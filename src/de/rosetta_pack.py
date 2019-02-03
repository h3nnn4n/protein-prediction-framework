import pyrosetta
import sys
import random
import time
import datetime
import string
import math
import bounds
import uuid
import protein_loader
import os

sys.path.append('../external')

import ramachandran.src.ramachandran as rama  # pylint: disable=E0401
import stride.src.stride as stride_  # pylint: disable=E0401
import tmscore.src.TMscore as tmscore  # pylint: disable=E0401


INIT = False

if not INIT:  # Initialize rosetta only once
    pyrosetta.init('-out:level 0 -ignore_unrecognized_res')

    INIT = True


"""
    no quaternary strut:
        - 1crn
        - 1zdd
        - 1enh

    quaternary strut:
        - 1rop
        - 1utg
        - 1ail
"""


class RosettaPack():
    def __init__(self, name='1zdd'):
        self.name = name

        self.allatom_switch = pyrosetta.SwitchResidueTypeSetMover('fa_standard')
        self.centroid_switch = pyrosetta.SwitchResidueTypeSetMover('centroid')
        self.scorefxn = pyrosetta.get_fa_scorefxn()

        self.scores = {}
        # self.scores['score0'] = pyrosetta.create_score_function('score0')
        # self.scores['score1'] = pyrosetta.create_score_function('score1')
        # self.scores['score2'] = pyrosetta.create_score_function('score2')
        # self.scores['score3'] = pyrosetta.create_score_function('score3')
        # self.scores['score5'] = pyrosetta.create_score_function('score5')
        # self.scores['scorefxn'] = pyrosetta.get_fa_scorefxn()

        self.ramachandran = rama.Ramachandran()
        p1 = "/usr/local/lib/python3.5/dist-packages/pyrosetta-4.0-py3.5-linux-x86_64.egg/pyrosetta/database/scoring/score_functions/" + \
             "rama/shapovalov/kappa75/all.ramaProb"
        p2 = "/usr/lib/python3.5/site-packages/pyrosetta-4.0-py3.5-linux-x86_64.egg/pyrosetta/database/scoring/score_functions/rama/s" + \
             "hapovalov/kappa75/all.ramaProb"
        p3 = "/usr/local/lib/python3.6/site-packages/pyrosetta-2018.22+release.99b36feae43-py3.6-macosx-10.13-x86_64.egg/pyrosetta/da" + \
             "tabase/scoring/score_functions/rama/shapovalov/kappa75/all.ramaProb"

        if os.path.exists(p1):
            self.ramachandran.load(path=p1)
        elif os.path.exists(p2):
            self.ramachandran.load(path=p2)
        else:
            self.ramachandran.load(path=p3)

        self.ramachandran.process()

        self.stride = stride_.Stride()
        self.stride.change_path("/home/h3nnn4n/Downloads/stride")

        self.tmscore = tmscore.TMscore(path='../external/tmscore/src/TMscore')

        self.fragset3 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(3)
        self.fragset9 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)

        self.bounds = bounds.Bounds()

        ######################

        self.protein_loader = protein_loader.ProteinLoader()
        self.protein_loader.load(name)
        self.name, self.target, self.ss_pred, self.native_path, self.fragset3_path, self.fragset9_path = self.protein_loader.get_data()

        self.native = pyrosetta.pose_from_pdb(self.native_path)

        self.fragset3.read_fragment_file(self.fragset3_path)
        self.fragset9.read_fragment_file(self.fragset9_path)

        ######################

        self.movemap = pyrosetta.MoveMap()
        self.movemap.set_bb(True)

        ######################

        self.pymover = pyrosetta.PyMOLMover()

        self.pose = pyrosetta.pose_from_sequence(self.target)
        self.loop_pose_cent = pyrosetta.pose_from_sequence(self.target)

        self.centroid_switch.apply(self.pose)

        self.pose.pdb_info().name(name)

        temp = 1.0
        n_moves = 2

        self.smallmover = pyrosetta.rosetta.protocols.simple_moves.SmallMover(self.movemap, temp, n_moves)
        self.shearmover = pyrosetta.rosetta.protocols.simple_moves.ShearMover(self.movemap, temp, n_moves)

        self.minmover = self.get_new_min_mover()

        cost = pyrosetta.rosetta.protocols.simple_moves.GunnCost()

        self.mover_3mer = pyrosetta.rosetta.protocols.simple_moves.ClassicFragmentMover(self.fragset3, self.movemap)
        self.mover_9mer = pyrosetta.rosetta.protocols.simple_moves.ClassicFragmentMover(self.fragset9, self.movemap)
        self.mover_3mer_smooth = pyrosetta.rosetta.protocols.simple_moves.SmoothFragmentMover(self.fragset3, self.movemap, cost)
        self.mover_9mer_smooth = pyrosetta.rosetta.protocols.simple_moves.SmoothFragmentMover(self.fragset9, self.movemap, cost)

        self.mover_3mer.set_movemap(self.movemap)
        self.mover_9mer.set_movemap(self.movemap)
        self.mover_3mer_smooth.set_movemap(self.movemap)
        self.mover_9mer_smooth.set_movemap(self.movemap)

        self.task_pack = pyrosetta.standard_packer_task(self.pose)
        self.task_pack.restrict_to_repacking()
        self.task_pack.or_include_current(True)
        self.pack_mover = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(self.scorefxn, self.task_pack)

        self.fast_relax = pyrosetta.rosetta.protocols.relax.FastRelax(self.scorefxn)

        self.minmover.movemap(self.movemap)
        self.minmover.score_function(self.get_score_function('score3'))

        self.mc = pyrosetta.MonteCarlo(self.pose, self.get_score_function('score0'), 2.0)

        self.abinitio = pyrosetta.rosetta.protocols.abinitio.ClassicAbinitio(self.fragset3, self.fragset9, self.movemap)
        self.abinitio.init(self.pose)

    def random_rama_angles_to_pose(self, pose):
        for k, v in enumerate(self.target):
            phi, psi = self.ramachandran.get_random_angle(v)
            pose.set_phi(k + 1, phi)
            pose.set_psi(k + 1, psi)

    def get_new_pose(self):
        pose = pyrosetta.pose_from_sequence(self.target)
        self.random_rama_angles_to_pose(pose)
        self.centroid_switch.apply(pose)
        return pose

    def get_rmsd_from_pose(self, pose=None, pose2=None):
        if pose2 is None:
            raise TypeError('Rosetta pack with internal states is no longer suported')
        else:
            return pyrosetta.rosetta.core.scoring.CA_rmsd(pose, pose2)

    def get_native(self):
        return self.native

    def get_rosetta_abinitio_protocol(self):
        return self.abinitio

    def run_tmscore(self, name=None, pose=None):
        clean_after = False
        if name is None and pose is not None:
            clean_after = True
            name = self.dump_tmp(pose)
        elif name is None and pose is None:
            raise ValueError('The name and pose arguments cant both be None')
        elif name is not None and pose is not None:
            raise ValueError('The name and pose arguments cant both be set')

        self.tmscore(self.native_path, name)

        if clean_after:
            os.remove(name)

    def get_tmscore(self):
        return self.tmscore.get_all()

    def get_sidechain_recover(self):
        return pyrosetta.rosetta.protocols.simple_moves.ReturnSidechainMover

    def get_min_mover(self):
        return self.minmover

    def get_new_min_mover(self, mm=None, score='score3_smooth',
                          alg='dfpmin_armijo_nonmonotone', tol=0.001):
        if mm is None:
            mm = self.movemap
        return pyrosetta.rosetta.protocols.minimization_packing.MinMover(
            mm,
            self.get_score_function(score),
            alg,
            tol,
            True
        )

    def get_small_mover(self):
        return self.smallmover

    def get_shear_mover(self):
        return self.shearmover

    def get_new_small_mover(self, temp=2.0, n_moves=1):
        return pyrosetta.rosetta.protocols.simple_moves.SmallMover(self.movemap, temp, n_moves)

    def get_new_shear_mover(self, temp=2.0, n_moves=1):
        return pyrosetta.rosetta.protocols.simple_moves.ShearMover(self.movemap, temp, n_moves)

    def get_9mer(self):
        return self.mover_9mer

    def get_3mer(self):
        return self.mover_3mer

    def get_9mer_smooth(self):
        return self.mover_9mer_smooth

    def get_3mer_smooth(self):
        return self.mover_3mer_smooth

    def get_mc(self):
        return self.mc

    def get_new_mc(self, pose, score, temp=2.0):
        return pyrosetta.MonteCarlo(pose, score, temp)

    def get_new_seq_mover(self):
        return pyrosetta.rosetta.protocols.moves.SequenceMover()

    def get_new_trial_mover(self, seq, mc):
        return pyrosetta.TrialMover(seq, mc)

    def get_new_rep_mover(self, trial, n):
        return pyrosetta.rosetta.protocols.moves.RepeatMover(trial, n)

    def get_packer(self):
        return self.pack_mover

    def get_fast_relax(self):
        return self.fast_relax

    def get_centroid_switch(self):
        return self.centroid_switch

    def get_allatom_switch(self):
        return self.allatom_switch

    def copy_pose_to_allatom(self, pose):
        ap = pyrosetta.Pose()
        ap.assign(pose)

        if ap.is_centroid():
            self.allatom_switch.apply(ap)

        return ap

    def convert_to_allatom_pose(self, pose):
        if pose.is_centroid():
            self.allatom_switch.apply(pose)

    def get_score_function(self, name='score3'):
        if name not in self.scores.keys():
            if name == 'scorefxn':
                self.scores[name] = pyrosetta.get_fa_scorefxn()
            else:
                self.scores[name] = pyrosetta.create_score_function(name)

        return self.scores[name]

    def get_score0(self):
        return self.get_score_function('score0')

    def get_score1(self):
        return self.get_score_function('score1')

    def get_score2(self):
        return self.get_score_function('score2')

    def get_score3(self):
        return self.get_score_function('score3')

    def get_score5(self):
        return self.get_score_function('score5')

    def get_scorefxn(self):
        return self.get_score_function('scorefxn')

    def get_new_movemap(self, free_backbone=True):
        movemap = pyrosetta.MoveMap()
        movemap.set_bb(free_backbone)
        return movemap

    def set_movemap_to_coil_and_loop_only(self, movemap):
        movemap.set_bb(False)
        for k, ss in enumerate(self.ss_pred):
            if ss == 'L' or ss == 'C':
                movemap.set_bb(k + 1, True)

    def get_new_movemap_with_free_coil_and_loop(self):
        movemap = self.get_new_movemap()
        self.set_movemap_to_coil_and_loop_only(movemap)
        return movemap

    def set_pose(self, pose, data):
        for i in range(len(self.target)):
            pose.set_phi(i + 1, data[i * 3 + 0])
            pose.set_psi(i + 1, data[i * 3 + 1])
            pose.set_omega(i + 1, data[i * 3 + 2])

    def get_backbone_angles_from_pose(self, pose):
        angles = []
        for i in range(len(self.target)):
            angles.append(pose.phi(i + 1))
            angles.append(pose.psi(i + 1))
            angles.append(pose.omega(i + 1))

        return angles

    def dump_tmp(self, pose):
        name = os.getcwd() + '/' + uuid.uuid4().hex
        pose.dump_pdb("%s" % (name))
        return name

    def dump(self, name, score, rmsd, prefix):
        now = datetime.datetime.now()
        char_set = string.ascii_uppercase + string.digits
        r_string = ''.join(random.sample(char_set * 6, 6))

        name_prefix = ""
        if prefix is not None:
            name_prefix = prefix + "_"

        name_sufix = "_%s_%04d_%02d_%02d__%02d_%02d_%02d__%s" % (name, now.year, now.month, now.day, now.hour, now.minute,
                                                                 now.second, r_string)

        self.pose.dump_pdb("%s%s_%010.5f_%010.5f_%s.pdb" %
                           (name_prefix, name, score, rmsd, name_sufix))
