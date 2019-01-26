import pytest
import os
from rosetta_pack import RosettaPack

pname = '1zdd'
eps = 1e-8

def test_path_doesnt_change():
    """
        Test that the working dir does not change
    """
    original_path = os.getcwd()
    RosettaPack(name=pname)
    after_path = os.getcwd()
    RosettaPack(name=pname)
    after_path2 = os.getcwd()

    assert original_path == after_path
    assert original_path == after_path2


def test_fragset_loaded():
    """
        Test that the 3 and 9mers fragments were loaded
    """

    rp = RosettaPack(name=pname)

    assert rp.fragset3.size() > 0
    assert rp.fragset9.size() > 0


def test_random_rama_angles_to_pose():
    """
        Test that random_rama_angles_to_pose changes all angles
    """

    rp = RosettaPack(name=pname)
    pose = rp.get_new_pose()

    old_angles = rp.get_backbone_angles_from_pose(pose)
    rp.random_rama_angles_to_pose(pose)
    new_angles = rp.get_backbone_angles_from_pose(pose)

    for k, (a, b) in enumerate(zip(new_angles, old_angles)):
        if (k + 1) % 3 == 0:
            continue  # Skip omega angles
        assert abs(a - b) > eps, "Failed at angle #%d" % k


def test_get_new_pose():
    rp = RosettaPack(name=pname)
    pose1 = rp.get_new_pose()
    pose2 = rp.get_new_pose()

    assert pose1 is not pose2


def test_get_rmsd_from_pose():
    rp = RosettaPack(name=pname)
    pose1 = rp.get_new_pose()
    pose2 = rp.get_new_pose()

    assert rp.get_rmsd_from_pose(pose1, pose2) > 1
    assert rp.get_rmsd_from_pose(pose1) > 1


def test_get_native():
    rp = RosettaPack(name=pname)
    pose = rp.get_new_pose()
    native = rp.get_native()

    assert pose is not native
    assert native.sequence() == rp.target


def test_run_tmscore():
    rp = RosettaPack(name=pname)
    rp.run_tmscore()
    data = rp.get_tmscore()

    assert len(data['gdt_ha']) == 2
    assert len(data['gdt_ha'][1]) == 4

    assert len(data['gdt_ts']) == 2
    assert len(data['gdt_ts'][1]) == 4

    assert data['tm_score'] is not None
    assert data['maxsub'] is not None
    assert data['rmsd'] is not None


def test_get_sidechain_recover():
    # TODO
    rp = RosettaPack(name=pname)
    rp.get_sidechain_recover()


def test_get_min_mover():
    rp = RosettaPack(name=pname)
    pose = rp.get_new_pose()
    min_mover = rp.get_min_mover()
    score_function = rp.get_score_function()
    rp.set_pose(pose, [0 for _ in range(3 * len(pose.sequence()))])

    score_before = score_function(pose)
    min_mover.apply(pose)
    score_after = score_function(pose)

    assert score_after < score_before


def get_new_min_mover():
    rp = RosettaPack(name=pname)
    pose = rp.get_new_pose()
    min_mover = rp.get_min_mover()
    min_mover2 = rp.get_new_min_mover()
    min_mover3 = rp.get_new_min_mover()
    score_function = rp.get_score_function()

    score_before = score_function(pose)
    min_mover2.apply(pose)
    score_after = score_function(pose)

    assert score_after < score_before

    assert min_mover is not min_mover2
    assert min_mover is not min_mover3
    assert min_mover2 is not min_mover3


def test_get_small_mover():
    rp = RosettaPack(name=pname)
    pose = rp.get_new_pose()
    mover = rp.get_small_mover()

    old_angles = rp.get_backbone_angles_from_pose(pose)
    mover.apply(pose)
    new_angles = rp.get_backbone_angles_from_pose(pose)

    assert old_angles != new_angles


def get_new_small_mover():
    rp = RosettaPack(name=pname)
    pose = rp.get_new_pose()
    small_mover = rp.get_small_mover()
    small_mover2 = rp.get_new_small_mover()
    small_mover3 = rp.get_new_small_mover()

    old_angles = rp.get_backbone_angles_from_pose(pose)
    small_mover2.apply(pose)
    new_angles = rp.get_backbone_angles_from_pose(pose)

    assert old_angles != new_angles

    assert small_mover is not small_mover2
    assert small_mover is not small_mover3
    assert small_mover2 is not small_mover3


def test_get_shear_mover():
    rp = RosettaPack(name=pname)
    pose = rp.get_new_pose()
    mover = rp.get_shear_mover()

    old_angles = rp.get_backbone_angles_from_pose(pose)
    mover.apply(pose)
    new_angles = rp.get_backbone_angles_from_pose(pose)

    assert old_angles != new_angles


def get_new_shear_mover():
    rp = RosettaPack(name=pname)
    pose = rp.get_new_pose()
    shear_mover = rp.get_shear_mover()
    shear_mover2 = rp.get_new_shear_mover()
    shear_mover3 = rp.get_new_shear_mover()

    old_angles = rp.get_backbone_angles_from_pose(pose)
    shear_mover2.apply(pose)
    new_angles = rp.get_backbone_angles_from_pose(pose)

    assert old_angles != new_angles

    assert shear_mover is not shear_mover2
    assert shear_mover is not shear_mover3
    assert shear_mover2 is not shear_mover3


def test_fragments():
    rp = RosettaPack(name=pname)
    pose = rp.get_new_pose()

    mover3mer = rp.get_3mer()
    mover9mer = rp.get_9mer()
    mover3mer_smooth = rp.get_3mer_smooth()
    mover9mer_smooth = rp.get_9mer_smooth()

    movers = [mover3mer, mover9mer, mover3mer_smooth, mover9mer_smooth]

    for mover in movers:
        old_angles = rp.get_backbone_angles_from_pose(pose)
        mover.apply(pose)
        new_angles = rp.get_backbone_angles_from_pose(pose)
        assert old_angles != new_angles

    for i in range(len(movers)):
        for j in range(i + 1, len(movers)):
            assert movers[i] is not movers[j]


def test_set_pose():
    rp = RosettaPack(name=pname)
    pose = rp.get_new_pose()

    angles = rp.get_backbone_angles_from_pose(pose)
    rp.random_rama_angles_to_pose(pose)
    angles2 = rp.get_backbone_angles_from_pose(pose)
    rp.set_pose(pose, angles)
    angles3 = rp.get_backbone_angles_from_pose(pose)

    assert angles != angles2
    assert angles == angles3
    assert angles2 != angles3
