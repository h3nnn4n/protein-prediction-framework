import pickle
import os

import numpy as np

    
def column_output(results):
    columns = ['rmsd_before', 'repack_time', 'score', 'rmsd_after', 'scorefxn']
    output = {}
    
    
    for column in columns:
        output[column] = []
        
        for result in results:
            output[column].append(result[column])

    return output


def experiment_summary(data, mode='best_by_energy', metric='scorefxn', with_raw=False):
    output = {}
    proteins = list(data.keys())
    
    for protein in proteins:
        output[protein] = {}

        experiments = list(data[protein].keys())
        
        for experiment in experiments:
            output[protein][experiment] = {}
            output[protein][experiment]['data'] = {}
            output[protein][experiment]['keys'] = ['min', 'max', 'mean', 'std', 'median']
            
            results = data[protein][experiment][mode]

            raw_data = column_output(results)[metric]
            
            if len(raw_data) == 0:
                print('[WARN] experiment %s for protein %s is empty! skipping' % (experiment, protein))
                continue
            
            output[protein][experiment]['data']['min'] = min(raw_data)
            output[protein][experiment]['data']['max'] = max(raw_data)
            output[protein][experiment]['data']['mean'] = np.mean(raw_data)
            output[protein][experiment]['data']['std'] = np.std(raw_data)
            output[protein][experiment]['data']['median'] = np.median(raw_data)
            
            if with_raw:
                output[protein][experiment]['data']['raw'] = raw_data

    return output

def load_all_data(runs):
    dataset = []
    
    for run in runs:
        os.chdir(run)
            
        experiment_folders = [file for file in os.listdir() if run in file]
        
        if len(experiment_folders) > 1:
            print('WARN: Found %d experiment folders for %s' % (len(experiment_folders), run))
        
        os.chdir(experiment_folders[-1])
        
        p_data = pickle.load(open( "data_dump.pickle", "rb"))
        dataset.append(p_data)
        os.chdir('../../')
        
    print('INFO: Loaded %d experiment runs dataset' % len(dataset))    
    
    return dataset


def merge_container(a, b):
    if type(a) == list:
        a.extend(b)
    elif type(a) == dict:
        a.update(b)
        
    return a
    
    
def merge_data(dataset, protein_blacklist=[], keep_only_common_methods=True):
    merged = {}
    
    modes = ['best_by_rmsd', 'all_repacks', 'best_by_energy']
    
    for data in dataset:
        for protein_k, protein_v in data.items():
            if protein_k not in merged.keys():
                merged[protein_k] = {}
                
            for experiment_k, experiment_v in protein_v.items():
                if experiment_k not in merged[protein_k].keys():
                    merged[protein_k][experiment_k] = {}
                    
                for mode in modes:
                    if mode not in merged[protein_k][experiment_k].keys():
                        if mode == 'all_repacks':
                            merged[protein_k][experiment_k][mode] = {}
                        else:
                            merged[protein_k][experiment_k][mode] = []
                        
                    merge_container(merged[protein_k][experiment_k][mode], experiment_v[mode])
    
    if len(protein_blacklist) > 0:
        remove_blacklisted_proteins(merged, protein_blacklist)
        
    if keep_only_common_methods:
        remove_not_common_methods(merged)
    
    return merged


def remove_blacklisted_proteins(alldata, protein_blacklist):
    proteins_before = list(alldata.keys())

    for protein in protein_blacklist:
        alldata.pop(protein.lower(), False)

    proteins_after = list(alldata.keys())

    print('removed %d proteins. Blacklist had %d' % (
        len(proteins_before) - len(proteins_after),
        len(protein_blacklist)
    ))


def remove_not_common_methods(alldata):
    first_key = list(alldata.keys())[0]
    methods = set(alldata[first_key].keys())

    for _, protein_methods in alldata.items():
        methods_ = set(protein_methods.keys())
        methods = methods & methods_

    proteins = list(alldata.keys())
    for protein in proteins:
        protein_methods = list(alldata[protein].keys())

        for protein_method in protein_methods:
            if protein_method not in methods:
                print('[WARN] removing %s from %s' % (protein_method, protein))
                alldata[protein].pop(protein_method, True)


def flatten_all_repacks_data(data):
    for protein_k, protein_v in data.items():
        for experiment_k, experiment_v in protein_v.items():
            source = experiment_v['all_repacks'].values()
            flattened = []
            
            for s in source:
                flattened.extend(s)
                
            experiment_v['all_repacks'] = flattened
    
    return data