import pytest
import mock
import os
from rosetta_classic_abinitio import ClassicAbinitio, main


def test_it_runs():
    files_before = os.listdir()
    ClassicAbinitio(pname='1enh').run(factor=1)
    files_after = os.listdir()

    s = set(files_before)
    new_files = [x for x in files_after if x not in s and '1enh_classic-abinitio_' in x]

    assert prefix_in_filelist('best_0', new_files)
    assert prefix_in_filelist('best_repacked_', new_files)
    assert prefix_in_filelist('stats_', new_files)
    assert prefix_in_filelist('repack_', new_files)
    assert prefix_in_filelist('parameters_', new_files)

    for file in new_files:
        assert os.path.getsize(file) > 0

    for file in new_files:
        os.remove(file)


def test_it_runs_as_script():
    """
        Mock dump_pdb_repacked so the test runs faster. It is (or should be) already tested elsewhere
    """
    mocked_sys_argv = ['rosetta_classic_abinitio.py']

    with mock.patch('rosetta_classic_abinitio.ClassicAbinitio.dump_pdb_repacked', mock.MagicMock(return_value=42)):
        with mock.patch('sys.argv', mocked_sys_argv):
            files_before = os.listdir()
            main()
            files_after = os.listdir()

            s = set(files_before)
            new_files = [x for x in files_after if x not in s and '1crn_classic-abinitio_' in x]

            assert prefix_in_filelist('best_0', new_files)
            # The following isnt generated because the function that logs it was mocked out
            # assert prefix_in_filelist('best_repacked_', new_files)
            assert prefix_in_filelist('stats_', new_files)
            # assert prefix_in_filelist('repack_', new_files)
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
