import datetime
import tarfile
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
    
    
def move_to_results(target):
    char_set = string.ascii_uppercase + string.digits
    r_string = ''.join(random.sample(char_set * 6, 6))

    now = datetime.datetime.now()

    name_suffix = "__%04d_%02d_%02d" % (
        now.year, now.month, now.day
    )
    
    folder_name = target + name_suffix
    
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
        print('INFO: created results folder %s' % folder_name)
    
    protein_list = get_protein_list()
    
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

    file_counter = 0
    total = len(os.listdir())
    file_counter = 0

    print('INFO: found %d files' % total)

    for filename in os.listdir():
        if '.yaml' in filename or '.dat' in filename or '.pdb': # in filename or 'spicker_' in filename:
            if any(map(lambda text: text in filename, BLACKLIST)):
                continue

            if any(map(lambda text: text in filename, WHITELIST)):
                file_counter += 1

                if file_counter % 1000 == 0:
                    print('INFO: Moved %7d of %7d' % (file_counter, total))

                protein = detect_protein(filename, protein_list)
                
                move(filename, os.path.join(folder_name, protein, filename))

    print('INFO: Moved %7d of %7d files' % (file_counter, total))

    return folder_name


def organize_protein_folder(protein):
    os.chdir(protein)
    
    print('INFO: Organizing %s' % protein)
    
    for target in WHITELIST:
        folder_name = target.replace('_', '')
        
        if not os.path.exists(folder_name):
            os.mkdir(folder_name)
        
        for file in os.listdir():
            if target in file:
                move(file, os.path.join(folder_name, file))
    
    os.chdir('..')

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
