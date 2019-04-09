import os
import mock
import pytest
from main import boot, main


mocked_config = {
    'pname': '1crn',
    'pop_size': 50,
    'max_evals': 1000,
    'stop_condition': 'evals',
    'stage0_init': 1,
    'sade_run': 1,
    'sade_lp': 10,
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


def test_run_from_main():
    files_before = os.listdir()
    with mock.patch('sys.argv', ['main.py']):
        main()
    files_after = os.listdir()

    s = set(files_before)
    new_files = [x for x in files_after if x not in s and '_none_' in x]

    assert prefix_in_filelist('best_0', new_files)
    assert prefix_in_filelist('best_repacked_', new_files)
    assert prefix_in_filelist('ops__', new_files)
    assert prefix_in_filelist('stats__', new_files)
    assert prefix_in_filelist('repack__', new_files)
    assert prefix_in_filelist('parameters__', new_files)
    assert prefix_in_filelist('forced_fragment__', new_files)

    for file in new_files:
        assert os.path.isfile(file)

    for file in new_files:
        os.remove(file)


def test_run_with_fake_file():
    with mock.patch("config_loader.ConfigLoader.get_config_file", mock.MagicMock(return_value=mocked_config)):
        files_before = os.listdir()
        boot('test-run')
        files_after = os.listdir()

        s = set(files_before)
        new_files = [x for x in files_after if x not in s and '_test-run_'in x]

        assert prefix_in_filelist('best_0', new_files)
        assert prefix_in_filelist('best_repacked_', new_files)
        assert prefix_in_filelist('ops__', new_files)
        assert prefix_in_filelist('stats__', new_files)
        assert prefix_in_filelist('repack__', new_files)
        assert prefix_in_filelist('parameters__', new_files)
        assert prefix_in_filelist('forced_fragment__', new_files)

        for file in new_files:
            assert os.path.isfile(file)

        for file in new_files:
            os.remove(file)


def test_run_lsh_with_fake_file():
    lsh_mocked_config = mocked_config.copy()
    lsh_mocked_config['do_lsh'] = True
    lsh_mocked_config['n_hashes'] = 2
    lsh_mocked_config['n_buckets'] = 5
    lsh_mocked_config['update_interval'] = 2
    lsh_mocked_config['change_interval'] = 5
    lsh_mocked_config['ops'] = ['rand1exp_lsh', 'monte_carlo_3']
    cname = 'test-lsh-run'

    with mock.patch("config_loader.ConfigLoader.get_config_file", mock.MagicMock(return_value=lsh_mocked_config)):
        files_before = os.listdir()
        boot(cname)
        files_after = os.listdir()

        s = set(files_before)
        new_files = [x for x in files_after if x not in s and ('_%s_' % cname) in x]

        assert prefix_in_filelist('best_0', new_files)
        assert prefix_in_filelist('best_repacked_', new_files)
        assert prefix_in_filelist('ops__', new_files)
        assert prefix_in_filelist('stats__', new_files)
        assert prefix_in_filelist('repack__', new_files)
        assert prefix_in_filelist('parameters__', new_files)
        assert prefix_in_filelist('forced_fragment__', new_files)

        for file in new_files:
            assert os.path.isfile(file)

        for file in new_files:
            os.remove(file)


# Utils

def prefix_in_filelist(prefix, filelist):
    for file in filelist:
        if prefix in file:
            return True

    return False
