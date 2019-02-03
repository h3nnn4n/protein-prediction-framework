import pytest
import os
import pyrosetta
from protein_loader import ProteinLoader

pname = '1crn'


def test_path_doesnt_change_on_init():
    """
        Test that the working dir does not change when __init__ is called
    """
    original_path = os.getcwd()
    ProteinLoader()
    after_path = os.getcwd()
    ProteinLoader()
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


def test_returns_correct_path():
    """
        Test that the relative paths to where the script was run are valid
        and that the primary structure (target) has the same length
        as the predicted secundary structure (ss_pred)
    """
    protein_loader = ProteinLoader()
    protein_loader.load(pname)
    name, target, ss_pred, native_path, fragset3_path, fragset9_path = protein_loader.get_data()

    assert name == pname
    assert os.path.isfile(native_path)
    assert os.path.isfile(fragset3_path)
    assert os.path.isfile(fragset9_path)

    assert len(target) == len(ss_pred)


def test_raises_for_file_not_found():
    """
        Test that FileNotFoundError is raised if the protein data is not found
    """
    protein_loader = ProteinLoader()

    with pytest.raises(FileNotFoundError):
        protein_loader.load('babana_file')


def test_that_it_runs():
    """
        Test that the protein loader gets the right files
    """
    pyrosetta.init('-out:level 0 -ignore_unrecognized_res')
    fragset3 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(3)
    fragset9 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)

    protein_loader = ProteinLoader()
    protein_loader.load(pname)
    _, target, _, native_path, fragset3_path, fragset9_path = protein_loader.get_data()

    native = pyrosetta.pose_from_pdb(native_path)

    fragset3.read_fragment_file(fragset3_path)
    fragset9.read_fragment_file(fragset9_path)

    assert native.sequence() == target
    assert fragset3.size() > 0
    assert fragset9.size() > 0


def test_that_it_runs_without_init_name():
    """
        Test that the protein loader gets the right files
    """
    pyrosetta.init('-out:level 0 -ignore_unrecognized_res')
    fragset3 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(3)
    fragset9 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)

    protein_loader = ProteinLoader(name=pname)
    protein_loader.load()
    _, target, _, native_path, fragset3_path, fragset9_path = protein_loader.get_data()

    native = pyrosetta.pose_from_pdb(native_path)

    fragset3.read_fragment_file(fragset3_path)
    fragset9.read_fragment_file(fragset9_path)

    assert native.sequence() == target
    assert fragset3.size() > 0
    assert fragset9.size() > 0


def test_load_raises_if_no_name_is_set():
    protein_loader = ProteinLoader()

    with pytest.raises(ValueError):
        protein_loader.load()
