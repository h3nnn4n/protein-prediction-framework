import os
import re
import sys

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def list_all_tests():
    old_path = os.getcwd()

    user_path = os.path.expanduser('~')
    base_path = os.path.join(user_path, 'progs/de_supimpa/src/de', 'tests_all')

    os.chdir(base_path)
    all_dirs = [file for file in os.listdir() if os.path.isdir(file)]

    os.chdir(old_path)

    return all_dirs


def set_experiment_path(experiment_path):
    user_path = os.path.expanduser('~')
    base_path = os.path.join(user_path, 'progs/de_supimpa/src/de', 'tests_all', experiment_path)

    os.chdir(base_path)
    print(os.getcwd())


def list_proteins():
    all_files = os.listdir()
    all_dirs = [file for file in all_files if os.path.isdir(file) and len(file) == 4]

    return all_dirs


def list_experiments():
    all_files = os.listdir()
    all_dirs = [file for file in all_files if os.path.isdir(file)]

    return all_dirs


def set_protein(pname):
    if not os.path.isdir(pname):
        os.chdir('..')

    if not os.path.isdir(pname):
        os.chdir('..')

    os.chdir(pname)

    print(os.getcwd())


def set_experiment(expname):
    all_dirs = [file for file in os.listdir() if os.path.isdir(file)]

    matches = [dir for dir in all_dirs if expname in dir]
    has_match = len(matches) > 0

    exact_match = os.path.isdir(expname)

    exact_match = os.path.isdir(expname)
    if not os.path.isdir(expname) and not has_match:
        print('going up')
        os.chdir('..')

    if exact_match:
        os.chdir(expname)
    elif has_match:
        os.chdir(matches[0])

    print(os.getcwd())


def list_stats():
    all_files = [file for file in os.listdir() if os.path.isfile(file) and 'stats__' in file]

    return all_files


def parse_num(num):
    try:
        return int(num)
    except BaseException:
        return float(num)


def read_data_from_file(filename):
    data = []
    with open(filename) as f:
        for line in f.readlines():
            tokens = re.sub(r"\s\s+", " ", line).strip().split(' ')
            numbers = list(map(lambda x: parse_num(x), tokens))

            data.append(numbers)

    return data


def get_metric_from_data_file(filename, index):
    data = []

    all_data = read_data_from_file(filename)

    for line in all_data:
        data.append(line[index])

    return data


def get_xy_from_data_file(filename, x, y):
    x_data = []
    y_data = []

    data = read_data_from_file(filename)

    for line in data:
        x_data.append(line[x])
        y_data.append(line[y])

    return x_data, y_data


def get_all_metric(index):
    all_data = []
    files = list_stats()

    for data_file in files:
        evals, data = get_xy_from_data_file(data_file, metrics['evals'], index)

        all_data.append({
            'evals': evals,
            'data': data
        })

    return all_data


def plot_all(metric):
    all_data = data = get_all_metric(metric)

    plt.figure(figsize=(16.5, 10))

    for exp_data in all_data:
        evals = exp_data['evals']
        data = exp_data['data']

        sns.lineplot(evals, data)

    # plt.ylim(0, 200)
    plt.xlim(0, 500000)


def average_all(metric):
    all_data = get_all_metric(metric)

    assert len(all_data) > 0, "all_data is [] !"

    target_evals = all_data[0]['evals']
    target_data = []

    assert len(target_evals) > 0, "target_evals is [] !"

    for index, target in enumerate(target_evals):
        target_data.append([])

        for data_pair in all_data:
            assert len(data_pair) > 0, "data_pair is [] !"

            evals, data = data_pair.values()
            evals, data = data_pair['evals'], data_pair['data']

            assert len(evals) > 0, "evals is [] !"
            assert len(data) > 0, "data is [] !"

            for i in range(1, len(evals)):
                a, b = evals[i - 1], evals[i]

                if a == target:
                    target_data[index].append(data[i - 1])
                    break

                if b == target:
                    target_data[index].append(data[i])
                    break

                if a < target and b > target:
                    target_data[index].append(lerp(data[i], data[i - 1], (target - a) / (b - a)))
                    break

            if index == len(target_evals) - 1:
                target_data[index].append(data[i])

    # print(list(map(len, target_data)))

    mean = [
        sum(x) / len(x) for x in target_data
    ]

    # print(len(target_evals), len(mean))

    l = min(len(target_evals), len(mean))

    target_evals = target_evals[0:l]
    mean = mean[0:l]

    return target_evals, mean


def plot_all_best_and_mean_score_and_avg_rmsd():
    evals, mean_score = average_all(metrics['mean_score'])
    _, best_score = average_all(metrics['best_score'])
    _, avg_rmsd = average_all(metrics['avg_rmsd'])
    # _, moment_of_inertia = average_all(metrics['moment_of_inertia'])
    # _, pedro_diver = average_all(metrics['pedro_diver'])

    fig, ax1 = plt.subplots(figsize=(16.5, 10))
    ax2 = ax1.twinx()

    ax1.plot(evals, mean_score, 'g-', label='mean_score')
    ax1.plot(evals, best_score, 'g', label='best_score', linestyle='-.')

    ax2.plot(evals, avg_rmsd, 'b-', label='avg_rmsd')

    ax1.set_xlabel('Evals')
    ax1.set_ylabel('Score3', color='g')
    ax2.set_ylabel('Diversity', color='b')

    ax1.set_xlim(0, 500000)
    ax1.set_ylim(bottom=0)
    ax2.set_ylim(bottom=0)

    # lines, labels = ax1.get_legend_handles_labels()
    # lines2, labels2 = ax2.get_legend_handles_labels()

    fig.legend(loc='upper right', bbox_to_anchor=(0.8, 0.8))

    plt.show()


def lerp(a, b, p):
    return (a * p) + ((1 - p) * b)


metrics = {
    'evals': 0,
    'iters': 1,
    'best_score': 2,
    'mean_score': 3,
    'pedro_diver': 4,
    'avg_rmsd': 5,
    'best_rmsd': 6,
    'moment_of_inertia': 7
}
