import pickle
import sys
import os
import numpy as np
from scikit_posthocs import posthoc_dunn


def load_data(target):
    original_path = None

    if os.path.isdir(target):
        original_path = os.getcwd()
        os.chdir(target)

    with open('%s.pickle' % target, 'rb') as f:
        data = pickle.load(f)

    if original_path is not None:
        os.chdir(original_path)

    return data


def dunn(data, alpha=0.05, metric=None):
    dunn_result = posthoc_dunn(list(data.values()), p_adjust='sidak')

    methods = list(data.keys())
    n_methods = len(methods)
    names = [method[5:][:-8] for method in methods]
    print("Comparing %s %s with %d runs each and alpha: %6.2f" % (metric, names, len(data[methods[0]]), alpha))

    for i in range(n_methods):
        method_a = methods[i]
        data_a = data[method_a]
        mean_a = np.mean(data_a)

        for j in range(n_methods):
            method_b = methods[j]
            data_b = data[method_b]
            mean_b = np.mean(data_b)

            p = dunn_result[i][j]
            if p < alpha / 2.0 and mean_a < mean_b:
                name_a = method_a[5:][:-8]
                name_b = method_b[5:][:-8]
                print("%25s  %25s  p: %8.4f  mean_a: %6.2f  mean_b: %6.2f" % (
                    name_a, name_b, p, mean_a, mean_b
                ))


def process(target):
    data = load_data(target)
    methods = list(data.keys())
    # metrics = ['rmsd_fxn', 'best_fxn', 'rmsd', 'best']
    metrics = ['rmsd_fxn', 'best_fxn']

    for metric in metrics:
        metric_data = {}
        for method in methods:
            metric_data[method] = data[method][metric]

        dunn(metric_data, 0.05, metric)
        print()


if __name__ == "__main__":
    targets = []

    if len(sys.argv) > 1:
        targets = sys.argv[1:]

    for target in targets:
        print(target)
        process(target)
