import pytest
import os
from protein_loader import ProteinLoader

pname = '1zdd'


def test_path_doesnt_change_on_init():
    """
        Test that the working dir does not change when __init__ is called
    """
    original_path = os.getcwd()
    pl = ProteinLoader()
    after_path = os.getcwd()
    pl = ProteinLoader()
    after_path2 = os.getcwd()

    assert original_path == after_path
    assert original_path == after_path2


def test_path_doesnt_change_on_load():
    """
        Test that the working dir does not change when load is called
    """
    original_path = os.getcwd()
    pl = ProteinLoader()
    pl.load(pname)
    after_path = os.getcwd()
    pl = ProteinLoader()
    pl.load(pname)
    after_path2 = os.getcwd()

    assert original_path == after_path
    assert original_path == after_path2
