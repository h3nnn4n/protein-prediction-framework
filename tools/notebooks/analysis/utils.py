import datetime
import tarfile
import pickle
import string
import random
import time
import sys
import os
import re

from shutil import copy, move


WHITELIST = [
    'stats_',
    'parameters_',
    'hooke-jeeves_',
    'forced_fragment_',
    'spicker__',
    'ops_',
    'best_',
    '_repacked_',
    '_repack_',
    'repack_'
]

BLACKLIST = [
    'base_'
]


def tar_dir(target):   
    tar_name = target + '__bkp.tar.gz'
    
    if os.path.exists(tar_name):
        print('WARN: %s already exists. Skipping' % tar_name)
        return
    
    print('INFO: tarballing to %s' % tar_name)
    
    start_time = time.time()
    
    tar = tarfile.open(tar_name, "w:gz")
    tar.add(target)
    tar.close()
    
    end_time = time.time()
    
    copy(tar_name, os.path.expanduser("~"))
    
    end_time2 = time.time()
    
    print('INFO: tarballing took %fs' % (end_time - start_time))
    print('INFO: tarball bkp took %fs' % (end_time2 - end_time))

    
def get_protein_list():
    protein_list = set()
    
    for file in os.listdir():
        if 'stats__' in file:
           protein_list.add(file[7:11])
        
        if len(file) == 4 and os.path.isdir(file):
            protein_list.add(file)
        
    return sorted(list(protein_list))


def detect_protein(name, protein_list):
    for protein in protein_list:
        if protein in name:
            return protein
        
    return None
    
    
def move_to_results(target, folder_name=None):
    if folder_name is None:
        char_set = string.ascii_uppercase + string.digits
        r_string = ''.join(random.sample(char_set * 6, 6))

        now = datetime.datetime.now()

        name_suffix = "__%04d_%02d_%02d" % (
            now.year, now.month, now.day
        )

        folder_name = target + name_suffix
    
    # Create protein folders
    
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
        print('INFO: created results folder %s' % folder_name)
    
    protein_list = get_protein_list()
    
    # Sort by protein
    
    for protein in protein_list:
        protein_path = os.path.join(folder_name, protein)
        if os.path.exists(protein_path):
            if os.path.isdir(protein_path):
                print('INFO: protein folder %s already exists' % protein_path)
            elif os.path.isfile(protein_path):
                print('WARN: protein folder %s is not a folder. wtf' % protein_path)
            else:
                print('WARN: protein folder %s has something funky' % protein_path)
        else:
            os.mkdir(protein_path)
            print('INFO: created protein folder for %s' % protein_path)

    # Sort by filetype
    
    experiments = get_experiment_names()
    print('INFO: Found %d experiments' % len(experiments))
    file_counter = 0
    total = len(os.listdir())
    file_counter = 0
    
    create_experiment_folders(folder_name, protein_list, experiments)

    print('INFO: found %d files' % total)

    for filename in os.listdir():
        if '.yaml' in filename or '.dat' in filename or '.pdb': # in filename or 'spicker_' in filename:
            if any(map(lambda text: text in filename, BLACKLIST)):
                continue

            if any(map(lambda text: text in filename, WHITELIST)):
                file_counter += 1

                if file_counter % 1000 == 0:
                    print('INFO: Moved %7d of %7d' % (file_counter, total))
                    
                experiment_name = get_experiment_name(experiments, filename)

                protein = detect_protein(filename, protein_list)
                
                move(filename, os.path.join(folder_name, protein, experiment_name, filename))

    print('INFO: Moved %7d of %7d files' % (file_counter, total))

    return folder_name


def get_experiment_name(experiment_names, file):
    matches = []
    
    for name in experiment_names:
        if name in file:
            matches.append(name)
            
    if len(matches) == 0:
        return None
    
    if len(matches) == 1:
        return matches[0]
    
    return sorted(matches, key=len)[-1]


def create_experiment_folders(folder_name, proteins, experiments):
    for protein in proteins:
        for experiment in experiments:
            path = os.path.join(folder_name, protein, experiment)

            if os.path.exists(path):
                print('INFO: Experiment %s for %s already exists. Skipping' % (experiment, protein))
            else:
                os.mkdir(path)


def get_experiment_names():
    experiments = set()
    
    for file in os.listdir():
        if 'stats__' not in file:
            continue
            
        experiment_name = file.split('__')[2][5:]
        experiments.add(experiment_name)
        
    return sorted(list(experiments))


def organize_protein_folder(protein):
    os.chdir(protein)
    
    print('INFO: Organizing %s' % protein)
    
    experiments = [file for file in os.listdir() if os.path.isdir(file) and 'spicker___' not in file]
    
    for experiment in experiments:
        os.chdir(experiment)
        
        for target in WHITELIST:
            folder_name = target.replace('_', '')

            if not os.path.exists(folder_name):
                os.mkdir(folder_name)

            for file in os.listdir():
                if target in file:
                    move(file, os.path.join(folder_name, file))
                    
        os.chdir('..')
    
    os.chdir('..')


    # More like some repack data
def extract_all_repack_data():
    data = {}
    wanted_data = ['repack_time', 'score', 'scorefxn', 'rmsd_before', 'rmsd_after', 'gdt_ts_after', 'tm_score_after']
    repacked = [file for file in os.listdir() if '_repack_' in file or 'repack_' in file[0:7]]
    
    print('INFO: Parsing %7d repack.dat files' % len(repacked))

    for name in repacked:
        run_code = name.split('__')[-1][0:6]

        if run_code not in data.keys():
            data[run_code] = []

        with open(name, 'r') as f:
            new_data = {}

            for line in f.readlines():
                tokens = re.sub(' {2,}', ' ', line.strip()).split(' ')

                has_data = any([wanted == tokens[0][:-1] for wanted in wanted_data])
                if has_data:
                    new_data[tokens[0][:-1]] = float(tokens[1])
        
        data[run_code].append(new_data)
                    
    return data
      
    
def get_by_best_metric(data, metric):
    data_by_best = []
    
    for _, results in data.items():
        best = results[0]
        
        for result in results:
            if result[metric] < best[metric]:
                best = result
                
        data_by_best.append(best)
        
    return data_by_best

    
def get_by_best_rmsd():
    data = extract_all_repack_data()
    return get_by_best_metric(data, 'rmsd_after')


def get_by_best_energy():
    data = extract_all_repack_data()
    return get_by_best_metric(data, 'scorefxn')
    
    
def data_dump():
    data = {}

    protein_list = [file for file in os.listdir() if len(file) == 4 and os.path.isdir(file)]

    for protein in protein_list:
        data[protein] = {}

        experiments = os.listdir(protein)

        for experiment in experiments:
            os.chdir(os.path.join(protein, experiment, 'repack'))
            
            data[protein][experiment] = {}
            data[protein][experiment]['best_by_rmsd'] = get_by_best_rmsd()
            data[protein][experiment]['best_by_energy'] = get_by_best_energy()
            data[protein][experiment]['all_repacks'] = extract_all_repack_data()
            
            os.chdir('../../../')
            
    with open('%s.pickle' % 'data_dump', 'wb') as f:
        pickle.dump(data, f)
            
    return data
