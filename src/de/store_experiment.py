import sys
import os
from shutil import copyfile, move


WHITELIST = [
    'stats_',
    'parameters_',
    'hooke-jeeves_',
    'forced_fragment_',
    'ops_',
    'best_',
    '_repacked_',
    '_repack_',
    'repack_']


def check_target_folders():
    if not os.path.isdir('bkp_all'):
        os.mkdir('bkp_all')
        print('bkp_all not found! Creating it!')

    if not os.path.isdir('tests_all'):
        os.mkdir('tests_all')
        print('tests_all not found! Creating it!')


def create_exp_folders(name):
    if os.path.isdir(os.path.join('tests_all', name)):
        print('WARNING! A folder with the name `%s` already exists!' % name)
    else:
        os.mkdir(os.path.join('tests_all', name))

    if os.path.isdir(os.path.join('bkp_all', name)):
        print('WARNING! A folder with the name `%s` already exists!' % name)
    else:
        os.mkdir(os.path.join('bkp_all', name))


def file_list():
    for filename in os.listdir():
        if '.yaml' in filename or '.dat' in filename or '.pdb' in filename:
            if 'base_' in filename:
                continue

            if any(map(lambda text: text in filename, WHITELIST)):
                yield filename

        continue


def copy_and_move_file(exp_name, filename):
    copyfile(filename, os.path.join('tests_all', exp_name, filename))
    move(filename, os.path.join('bkp_all', exp_name, filename))


def main():
    args = sys.argv
    exp_name = args[1]

    print('creating experiment folder: %s' % exp_name)

    check_target_folders()
    create_exp_folders(exp_name)

    for filename in file_list():
        copy_and_move_file(exp_name, filename)


if __name__ == "__main__":
    main()
