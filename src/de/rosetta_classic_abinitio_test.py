import pytest
import os
from rosetta_classic_abinitio import ClassicAbinitio


def test_it_runs():
    files_before = os.listdir()
    ClassicAbinitio(pname='1crn').run(factor=1)
    files_after = os.listdir()

    s = set(files_before)
    new_files = [x for x in files_after if x not in s and 'classic-abinitio' in x]

    assert prefix_in_filelist('best_0', new_files)
    assert prefix_in_filelist('best_repacked_', new_files)
    assert prefix_in_filelist('stats_', new_files)
    assert prefix_in_filelist('repack_', new_files)
    assert prefix_in_filelist('parameters_', new_files)

    for file in new_files:
        assert os.path.getsize(file) > 0

    for file in new_files:
        os.remove(file)


# Utils

def prefix_in_filelist(prefix, filelist):
    for file in filelist:
        if prefix in file:
            return True

    return False
