import pytest
import os
from rosetta_pack import RosettaPack


pname = '1crn'
rp = RosettaPack(name=pname)
eps = 1e-8


def test_it_loads_more_than_once():
    rp = RosettaPack(name=pname)
    assert rp is not None
    rp = RosettaPack(name=pname)
    assert rp is not None

    import rosetta_pack

    rp = rosetta_pack.RosettaPack(name=pname)
    assert rp is not None
    rp = rosetta_pack.RosettaPack(name=pname)
    assert rp is not None

    import rosetta_pack

    rp = rosetta_pack.RosettaPack(name=pname)
    assert rp is not None
    rp = rosetta_pack.RosettaPack(name=pname)
    assert rp is not None


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
    assert rp.fragset3.size() > 0
    assert rp.fragset9.size() > 0


def test_random_rama_angles_to_pose():
    """
        Test that random_rama_angles_to_pose changes all angles
    """
    pose = rp.get_new_pose()

    old_angles = rp.get_backbone_angles_from_pose(pose)
    rp.random_rama_angles_to_pose(pose)
    new_angles = rp.get_backbone_angles_from_pose(pose)

    for k, (a, b) in enumerate(zip(new_angles, old_angles)):
        if (k + 1) % 3 == 0:
            continue  # Skip omega angles
        assert abs(a - b) > eps, "Failed at angle #%d" % k


def test_get_new_pose():
    pose1 = rp.get_new_pose()
    pose2 = rp.get_new_pose()

    assert pose1 is not pose2


def test_get_rmsd_from_pose():
    pose1 = rp.get_new_pose()
    pose2 = rp.get_new_pose()

    assert rp.get_rmsd_from_pose(pose1, pose2) > 1


def test_get_rmsd_from_pose_deprecation_raise():
    pose = rp.get_new_pose()

    with pytest.raises(TypeError):
        rp.get_rmsd_from_pose(pose)


def test_get_native():
    pose = rp.get_new_pose()
    native = rp.get_native()

    assert pose is not native
    assert native.sequence() == rp.target


def test_run_tmscore_with_pose():
    pose = rp.get_new_pose()
    rp.run_tmscore(pose=pose)
    data = rp.get_tmscore()

    assert len(data['gdt_ha']) == 2
    assert len(data['gdt_ha'][1]) == 4

    assert len(data['gdt_ts']) == 2
    assert len(data['gdt_ts'][1]) == 4

    assert data['tm_score'] is not None
    assert data['maxsub'] is not None
    assert data['rmsd'] is not None


def test_run_tmscore_with_path():
    pose = rp.get_new_pose()
    path = rp.dump_tmp(pose)
    rp.run_tmscore(name=path)
    data = rp.get_tmscore()

    assert len(data['gdt_ha']) == 2
    assert len(data['gdt_ha'][1]) == 4

    assert len(data['gdt_ts']) == 2
    assert len(data['gdt_ts'][1]) == 4

    assert data['tm_score'] is not None
    assert data['maxsub'] is not None
    assert data['rmsd'] is not None

    if os.path.isfile(path):
        os.remove(path)


def test_run_tmscore_with_no_args():
    with pytest.raises(ValueError):
        rp.run_tmscore()


def test_run_tmscore_with_too_many_args():
    pose = rp.get_new_pose()
    with pytest.raises(ValueError):
        rp.run_tmscore(name='potato', pose=pose)


def test_run_tmscore_path():
    original_path = os.getcwd()
    pose = rp.get_new_pose()
    rp.run_tmscore(pose=pose)
    after_path = os.getcwd()
    rp.run_tmscore(pose=pose)
    after_path2 = os.getcwd()

    assert original_path == after_path
    assert original_path == after_path2


def test_get_tmscore_path():
    original_path = os.getcwd()
    pose = rp.get_new_pose()
    rp.run_tmscore(pose=pose)
    rp.get_tmscore()
    after_path = os.getcwd()
    pose = rp.get_new_pose()
    rp.run_tmscore(pose=pose)
    rp.get_tmscore()
    after_path2 = os.getcwd()

    assert original_path == after_path
    assert original_path == after_path2


def test_get_sidechain_recover():
    # TODO
    rp = RosettaPack(name=pname)
    rp.get_sidechain_recover()


def test_get_min_mover():
    pose = rp.get_new_pose()
    min_mover = rp.get_min_mover()
    score_function = rp.get_score_function()
    rp.set_pose(pose, [0 for _ in range(3 * len(pose.sequence()))])

    score_before = score_function(pose)
    min_mover.apply(pose)
    score_after = score_function(pose)

    assert score_after < score_before


def test_get_new_min_mover():
    pose = rp.get_new_pose()
    min_mover = rp.get_min_mover()
    min_mover2 = rp.get_new_min_mover()
    min_mover3 = rp.get_new_min_mover(rp.movemap)
    score_function = rp.get_score_function()
    rp.set_pose(pose, [0 for _ in range(3 * len(pose.sequence()))])

    score_before = score_function(pose)
    min_mover2.apply(pose)
    score_after = score_function(pose)

    assert score_after < score_before

    assert min_mover is not min_mover2
    assert min_mover is not min_mover3
    assert min_mover2 is not min_mover3


def test_get_small_mover():
    pose = rp.get_new_pose()
    mover = rp.get_small_mover()

    old_angles = rp.get_backbone_angles_from_pose(pose)
    mover.apply(pose)
    new_angles = rp.get_backbone_angles_from_pose(pose)

    assert old_angles != new_angles


def test_get_new_small_mover():
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
    pose = rp.get_new_pose()
    mover = rp.get_shear_mover()

    old_angles = rp.get_backbone_angles_from_pose(pose)
    mover.apply(pose)
    new_angles = rp.get_backbone_angles_from_pose(pose)

    assert old_angles != new_angles


def test_get_new_shear_mover():
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
    pose = rp.get_new_pose()

    angles = rp.get_backbone_angles_from_pose(pose)
    rp.random_rama_angles_to_pose(pose)
    angles2 = rp.get_backbone_angles_from_pose(pose)
    rp.set_pose(pose, angles)
    angles3 = rp.get_backbone_angles_from_pose(pose)

    assert angles != angles2
    assert angles == angles3
    assert angles2 != angles3


def test_mc():
    pose = rp.get_new_pose()
    mc = rp.get_mc()
    mover9mer = rp.get_9mer()
    score = rp.get_score3()

    rp.set_pose(pose, [0 for _ in range(3 * len(pose.sequence()))])

    score_before = score(pose)

    for _ in range(10):
        mover9mer.apply(pose)

    score_after = score(pose)

    assert score_after < score_before
    assert mc.last_accepted_score() < score_before


def test_get_new_mc():
    pose = rp.get_new_pose()
    score = rp.get_score3()
    mc = rp.get_mc()
    mc2 = rp.get_new_mc(pose, score)
    mc3 = rp.get_new_mc(pose, score)
    mover9mer = rp.get_9mer()

    rp.set_pose(pose, [0 for _ in range(3 * len(pose.sequence()))])

    score_before = score(pose)

    for _ in range(10):
        mover9mer.apply(pose)

    score_after = score(pose)

    assert score_after < score_before
    assert mc.last_accepted_score() < score_before

    assert mc is not mc2
    assert mc is not mc3
    assert mc2 is not mc3


def test_get_new_trial_mover():
    pose = rp.get_new_pose()
    score = rp.get_score3()
    mc = rp.get_mc()
    mover9mer = rp.get_9mer()
    trial_mover = rp.get_new_trial_mover(mover9mer, mc)

    rp.set_pose(pose, [0 for _ in range(3 * len(pose.sequence()))])

    score_before = score(pose)

    for _ in range(10):
        trial_mover.apply(pose)

    score_after = score(pose)

    assert score_after < score_before


def test_get_new_seq_mover():
    pose = rp.get_new_pose()
    score = rp.get_score3()
    mc = rp.get_mc()
    mover3mer = rp.get_9mer()
    mover9mer = rp.get_9mer()

    seq_mover = rp.get_new_seq_mover()
    seq_mover.add_mover(mover3mer)
    seq_mover.add_mover(mover9mer)
    trial_mover = rp.get_new_trial_mover(seq_mover, mc)

    rp.set_pose(pose, [0 for _ in range(3 * len(pose.sequence()))])

    score_before = score(pose)

    for _ in range(10):
        trial_mover.apply(pose)

    score_after = score(pose)

    assert score_after < score_before


def test_get_new_rep_mover():
    pose = rp.get_new_pose()
    score = rp.get_score3()
    mc = rp.get_mc()
    mover9mer = rp.get_9mer()
    trial_mover = rp.get_new_trial_mover(mover9mer, mc)
    rep_mover = rp.get_new_rep_mover(trial_mover, 10)

    rp.set_pose(pose, [0 for _ in range(3 * len(pose.sequence()))])

    score_before = score(pose)
    rep_mover.apply(pose)
    score_after = score(pose)

    assert score_after < score_before


def test_get_centroid_switch():
    pose = rp.get_new_pose()
    switch = rp.get_centroid_switch()
    switch_allatom = rp.get_allatom_switch()

    assert pose.is_centroid()
    switch.apply(pose)
    assert pose.is_centroid()
    switch_allatom.apply(pose)
    assert pose.is_fullatom()
    switch.apply(pose)
    assert pose.is_centroid()


def test_get_allatom_switch():
    pose = rp.get_new_pose()
    switch = rp.get_allatom_switch()

    assert pose.is_centroid()
    switch.apply(pose)
    assert pose.is_fullatom()


def test_convert_to_allatom_pose():
    pose = rp.get_new_pose()
    assert pose.is_centroid()
    rp.convert_to_allatom_pose(pose)
    assert pose.is_fullatom()


def test_convert_to_allatom_pose_when_allatom():
    pose = rp.get_new_pose()
    rp.convert_to_allatom_pose(pose)
    assert pose.is_fullatom()
    rp.convert_to_allatom_pose(pose)
    assert pose.is_fullatom()


def test_copy_pose_to_allatom():
    pose = rp.get_new_pose()
    assert pose.is_centroid()
    new_pose = rp.copy_pose_to_allatom(pose)
    assert new_pose.is_fullatom()
    assert pose.is_centroid()
    assert pose is not new_pose


def test_copy_pose_to_allatom_when_allatom():
    pose = rp.get_new_pose()
    rp.convert_to_allatom_pose(pose)
    assert pose.is_fullatom()
    new_pose = rp.copy_pose_to_allatom(pose)
    assert new_pose.is_fullatom()
    assert pose is not new_pose


def test_get_new_movemap():
    movemap = rp.get_new_movemap()
    movemap2 = rp.get_new_movemap()

    assert movemap is not movemap2
    assert movemap is not rp.movemap
    assert movemap2 is not rp.movemap


def test_get_new_movemap_free_backbone():
    movemap = rp.get_new_movemap(free_backbone=True)
    movemap2 = rp.get_new_movemap(free_backbone=True)
    movemap3 = rp.get_new_movemap(free_backbone=False)

    for i in range(len(rp.target)):
        assert movemap.get_bb(i + 1)
        assert movemap2.get_bb(i + 1)
        assert not movemap3.get_bb(i + 1)


def test_set_movemap_to_coil_and_loop_only():
    movemap = rp.get_new_movemap()
    rp.set_movemap_to_coil_and_loop_only(movemap)

    for k, ss in enumerate(rp.ss_pred):
        if ss == 'L' or ss == 'C':
            assert movemap.get_bb(k + 1)
        else:
            assert not movemap.get_bb(k + 1)


def test_get_new_movemap_with_free_coil_and_loop():
    movemap = rp.get_new_movemap_with_free_coil_and_loop()

    for k, ss in enumerate(rp.ss_pred):
        if ss == 'L' or ss == 'C':
            assert movemap.get_bb(k + 1)
        else:
            assert not movemap.get_bb(k + 1)


def test_get_score0():
    scoref = rp.get_score0()
    scoref2 = rp.get_score0()

    assert scoref is scoref2
    assert scoref.get_name() == 'score0'


def test_get_score1():
    scoref = rp.get_score1()
    scoref2 = rp.get_score1()

    assert scoref is scoref2
    assert scoref.get_name() == 'score1'


def test_get_score2():
    scoref = rp.get_score2()
    scoref2 = rp.get_score2()

    assert scoref is scoref2
    assert scoref.get_name() == 'score2'


def test_get_score3():
    scoref = rp.get_score3()
    scoref2 = rp.get_score3()

    assert scoref is scoref2
    assert scoref.get_name() == 'score3'


def test_get_score5():
    scoref = rp.get_score5()
    scoref2 = rp.get_score5()

    assert scoref is scoref2
    assert scoref.get_name() == 'score5'


def test_get_scorefxn():
    scoref = rp.get_scorefxn()
    scoref2 = rp.get_scorefxn()

    assert scoref is scoref2
    assert scoref.get_name() == 'ref2015'


def test_get_rosetta_abinitio_protocol():
    ab1 = rp.get_rosetta_abinitio_protocol()
    ab2 = rp.get_rosetta_abinitio_protocol()

    assert ab1 is ab2

    pose = rp.get_new_pose()
    scoref = rp.get_score3()
    score_before = scoref(pose)
    ab1.apply(pose)
    score_after = scoref(pose)

    assert score_after < score_before


def test_get_fast_relax():
    fast_relax = rp.get_fast_relax()
    assert fast_relax.get_name() == 'FastRelax'


def test_get_packer():
    fast_relax = rp.get_packer()
    assert fast_relax.get_name() == 'PackRotamersMover'
