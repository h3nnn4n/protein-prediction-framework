import pytest
import os
from rosetta_pack import RosettaPack

pname = '1zdd'


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


def test_native_pose_loaded():
    """
        Test that pose was loaded with the native sequence
    """

    rp = RosettaPack(name=pname)

    assert rp.native.sequence() == rp.target
