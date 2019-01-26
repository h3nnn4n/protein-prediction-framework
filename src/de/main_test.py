import pytest
import os
from main import boot


def test_run():
    files_before = os.listdir()
    boot(None)
    files_after = os.listdir()

    s = set(files_before)
    new_files = [x for x in files_after if x not in s]

    assert prefix_in_filelist('best_0', new_files)
    assert prefix_in_filelist('best_repacked_', new_files)
    assert prefix_in_filelist('ops__', new_files)
    assert prefix_in_filelist('stats__', new_files)
    assert prefix_in_filelist('repack__', new_files)
    assert prefix_in_filelist('parameters__', new_files)

    for file in new_files:
        assert os.path.getsize(file) > 0

    for file in new_files:
        os.remove(file)


def prefix_in_filelist(prefix, filelist):
    for file in filelist:
        if prefix in file:
            return True

    return False
