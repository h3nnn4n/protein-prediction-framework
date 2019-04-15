import os
import re
import sys

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


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


def set_experiment(expname):
    if not os.path.isdir(expname):
        os.chdir('..')

    os.chdir(expname)


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

    # print(len(evals), len(mean))

    return evals, mean


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
