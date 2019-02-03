import pytest
import random
from protein_data import ProteinData
from rosetta_pack import RosettaPack

# TODO:
# Test side chain being putback after allatom -> centroid -> allatom
# Test TM-Score

pname = '1crn'
rp = RosettaPack(name=pname)
eps = 1e-8  # The tolerance for comparing floats
repeats = 10  # Number of times to repeat random tests
stage2_modes = ['3', '3s', '9', '9s', 'shear', 'small']


def test_init():
    """
        Just test that the class is correctly instantiated and
        wont throw any yikes
    """
    pd = ProteinData(rp, allatom=True)

    assert pd is not None  # Just asserting something


def test_angles_match_pose_after_init():
    """
        test if the angle vector matches the pose angles after init
    """
    for _ in range(repeats):
        pd = ProteinData(rp, allatom=True)

        for k, (a, b) in enumerate(get_angles_and_pose(pd)):
            assert abs(a - b) < eps, "failed at angle #%d" % k


def test_angles_match_pose_after_reset():
    """
        test if the angle vector matches the pose angles after reset
    """
    pd = ProteinData(rp, allatom=True)

    for _ in range(repeats):
        pd.reset()

        for k, (a, b) in enumerate(get_angles_and_pose(pd)):
            assert abs(a - b) < eps, "failed at angle #%d" % k


def test_angles_match_pose_after_new_angles():
    """
        test if the angle vector matches the pose angles after a new set of angles
    """
    pd = ProteinData(rp, allatom=True)

    for _ in range(repeats):
        angles = get_uniform_angles(size=pd.nsca)
        pd.new_angles(angles)

        for k, (a, b) in enumerate(get_angles_and_pose(pd)):
            assert abs(a - b) < eps, "failed at angle #%d" % k


def test_eval_scorefxn():
    """
        test if the score is working correctly
        This test here could use some improvements
    """
    pd = ProteinData(rp, allatom=True)

    for _ in range(repeats):
        pd.reset()
        angles_rama = pd.angles.copy()
        angles_random = get_uniform_angles(size=pd.nsca)
        angles_bad = get_flat_angles(size=pd.nsca)

        pd.new_angles(angles_rama)
        pd.eval()
        score_rama = pd.score

        pd.new_angles(angles_random)
        pd.eval()
        score_random = pd.score

        pd.new_angles(angles_bad)
        pd.eval()
        score_bad = pd.score

        assert score_rama < score_bad
        assert score_random < score_bad


def test_eval_score3_vs_scorefxn():
    """
        test if the scorefxn is bigger than score3, which it should pretty
        much always be
    """
    pd_centroid = ProteinData(rp, allatom=False)
    pd_all_atom = ProteinData(rp, allatom=True)

    for _ in range(repeats):
        pd_centroid.reset()
        pd_all_atom.reset()

        assert pd_centroid.score < pd_all_atom.score


def test_update_angle_from_pose():
    """
        test that the angles from pose are correctly copied to the angle vector
    """
    pd = ProteinData(rp, allatom=True)

    for _ in range(repeats):
        random_angles = get_uniform_angles(size=pd.nsca)
        set_pose_angles(pd, random_angles)
        pd.update_angle_from_pose()

        for k, (a, b) in enumerate(get_angles_and_pose(pd)):
            assert abs(a - b) < eps, "failed at angle #%d" % k


def test_score_improves_with_stage1_mc():
    """
        test that the score improves with stage1 mc
    """
    pd = ProteinData(rp, allatom=True)

    for _ in range(repeats):
        pd.new_angles(get_flat_angles(pd.nsca))
        pd.eval()
        score_before = pd.score

        pd.stage1_mc(n=100, temp=2.0)
        score_after = pd.score

        assert score_after < score_before, "Score should improve"


def test_angles_match_pose_after_stage1_mc():
    """
        test that angles match after stage1 mc
    """
    pd = ProteinData(rp, allatom=True)

    for _ in range(repeats):
        pd.new_angles(get_flat_angles(pd.nsca))
        pd.stage1_mc(n=100, temp=2.0)

        for k, (a, b) in enumerate(get_angles_and_pose(pd)):
            assert abs(a - b) < eps, "failed at angle #%d" % k


def test_stage2_mc_runs():
    """
        test that stage2 mc runs
    """
    pd = ProteinData(rp, allatom=True)

    for _ in range(repeats):
        for mode in stage2_modes:

            pd.new_angles(get_flat_angles(pd.nsca))
            pd.eval()

            pd.stage2_mc(n=10, temp=2.0, mode=mode)


def test_angles_match_pose_after_stage2_mc():
    """
        test that angles match after stage2 mc
    """
    pd = ProteinData(rp, allatom=True)

    for _ in range(repeats):
        for mode in stage2_modes:
            pd.new_angles(get_flat_angles(pd.nsca))
            pd.stage2_mc(n=10, temp=2.0)

            for k, (a, b) in enumerate(get_angles_and_pose(pd)):
                assert abs(a - b) < eps, "failed at angle #%d, mode: %s" % (k, mode)


def test_dont_raise_stage2_mc_when_wrong_function():
    """
        test that angles match after stage2 mc
    """
    pd = ProteinData(rp, allatom=True)
    pd.set_score_function('score3')
    pd.stage2_mc(n=10, temp=2.0)

    assert pd is not None


def test_angles_match_pose_after_copy():
    """
        test that angles match after copying another protein_data object
    """
    pd_origin = ProteinData(rp, allatom=True)
    pd_target = ProteinData(rp, allatom=True)

    for _ in range(repeats):
        pd_target.new_angles(get_flat_angles(pd_origin.nsca))
        pd_origin.reset()
        pd_target.copy(pd_origin)

        assert pd_origin.score == pd_target.score, "The score of both proteins should match"

        for k, (a, b) in enumerate(get_angles_and_pose(pd_origin)):
            assert abs(a - b) < eps, "failed at angle #%d" % k

        for k, (a, b) in enumerate(get_angles_and_pose(pd_target)):
            assert abs(a - b) < eps, "failed at angle #%d" % k

        for i in range(pd_origin.nsca):
            assert abs(pd_origin.angles[i] - pd_target.angles[i]) < eps, "Angle vector of both proteins should match"


def test_call_scores_make_sense():
    """
        simple test to verify if the value make some sense
    """
    pd = ProteinData(rp, allatom=True)

    for _ in range(repeats):
        pd.reset()
        angles_rama = pd.angles.copy()
        angles_random = get_uniform_angles(size=pd.nsca)
        angles_bad = get_flat_angles(size=pd.nsca)

        score_rama = pd(angles_rama)
        score_random = pd(angles_random)
        score_bad = pd(angles_bad)

        assert score_rama < score_bad
        assert score_random < score_bad


def test_call_evaluate_passed_angles():
    """
        test that the angles passed to __call__ are used to evaluate
    """
    pd = ProteinData(rp, allatom=True)

    for _ in range(repeats):
        pd.reset()
        angles_rama = pd.angles.copy()
        angles_random = get_uniform_angles(size=pd.nsca)

        pd.reset()  # Make sure that the rama angles dont match the ones in pose

        pd(angles_rama)
        for k, (a, b) in enumerate(zip(pd.angles, angles_rama)):
            assert abs(a - b) < eps, "failed at angle #%d" % k

        pd(angles_random)
        skip_count = 0
        for k, (a, b) in enumerate(zip(pd.angles, angles_random)):
            if a == 0:
                # Skip the first phi
                # and last psi and omega angles that are always zero
                skip_count += 1
                continue

            assert abs(a - b) < eps, "failed at angle #%d" % k

        assert skip_count == 3, "Should skip exactly 3 angles"


def test_repack():
    pd = ProteinData(rp)
    assert pd.repacked is None
    score = pd.repack()
    assert pd.repacked is not None
    assert type(score) is float

# Utils


def get_angles_and_pose(pd):
    angles = []
    pose_angles = []

    index = 0
    for k, _ in enumerate(pd.rosetta_pack.ss_pred):
        n_sidechain_angles = pd.bounds.getNumSideChainAngles(pd.rosetta_pack.target[k])

        phi_pose = pd.pose.phi(k + 1)
        psi_pose = pd.pose.psi(k + 1)
        omega_pose = pd.pose.omega(k + 1)

        pose_angles.append(phi_pose)
        pose_angles.append(psi_pose)
        pose_angles.append(omega_pose)

        phi_angle = pd.angles[index + 0]
        psi_angle = pd.angles[index + 1]
        omega_angle = pd.angles[index + 2]

        angles.append(phi_angle)
        angles.append(psi_angle)
        angles.append(omega_angle)

        if pd.allatom and n_sidechain_angles > 0:
                aa = pd.rosetta_pack.target[k]
                side_chain_angles = pd.bounds.generateRandomSidechainAngles(aa)
                for _, vv in enumerate(side_chain_angles):
                    angles.append(vv)
                    pose_angles.append(vv)

        index += 3 + n_sidechain_angles

    return zip(angles, pose_angles)


def get_uniform_angles(size=None):
    return [random.uniform(-180, 180) for _ in range(size)]


def get_flat_angles(size=None):
    return [0 for _ in range(size)]


def set_pose_angles(pd, angles):
    """
        Writes angles to the pose object only
    """
    index = 0
    size = len(pd.rosetta_pack.ss_pred)
    for k, _ in enumerate(pd.rosetta_pack.ss_pred):
        n_sidechain_angles = pd.bounds.getNumSideChainAngles(pd.rosetta_pack.target[k])
        phi, psi, omega = angles[index + 0], angles[index + 1], angles[index + 2]

        if k == 0:
            phi = 0.0
        elif k == size - 1:
            psi = 0.0
            omega = 0.0

        pd.pose.set_phi(k + 1, phi)
        pd.pose.set_psi(k + 1, psi)
        pd.pose.set_omega(k + 1, omega)

        if pd.allatom and n_sidechain_angles > 0:
            for j in range(n_sidechain_angles):
                vv = angles[index + 3 + j]
                pd.pose.set_chi(j + 1, k + 1, vv)

        index += 3 + n_sidechain_angles
