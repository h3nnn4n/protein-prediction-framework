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
