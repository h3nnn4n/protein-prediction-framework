import os
import mock
import pytest
from main import boot


mocked_config = {
    'pname': '1crn',
    'pop_size': 10,
    'max_evals': 500,
    'stop_condition': 'evals',
    'stage0_init': 1,
    'sade_run': 1,
    'sade_lp': 50,
    'sade_selection': 'roulette',
    'enable_remc': 0,
    'ops': ['rand1bin_global',
            'rand2bin_global',
            'best1bin_global',
            'best2bin_global',
            'currToBest_global',
            'currToRand_global',
            'currToRand_exp_global',
            'currToBest_exp_global',
            'rand1exp_global',
            'rand2exp_global',
            'best1exp_global',
            'best2exp_global',
            'monte_carlo_3',
            'monte_carlo_3s',
            'monte_carlo_9',
            'monte_carlo_9s',
            'monte_carlo_3_coil_only'],
    'energy_function': 'score3',
    'extended_diversity_measurements': False
}


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


def test_run_with_fake_file():
    with mock.patch("config_loader.ConfigLoader.get_config_file", mock.MagicMock(return_value=mocked_config)):
        boot('not_a_file')


# Utils

def prefix_in_filelist(prefix, filelist):
    for file in filelist:
        if prefix in file:
            return True

    return False
