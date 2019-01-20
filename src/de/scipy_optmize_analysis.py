import pickle
import numpy as np
from scikit_posthocs import posthoc_dunn


def load_data():
    with open('scipy_min_test.pickle', 'rb') as f:
        data = pickle.load(f)

    with open('scipy_min_dual_annealing.pickle', 'rb') as f:
        data = {**data, **pickle.load(f)}

    return data


data = load_data()

dunn_result = posthoc_dunn(list(data.values()), p_adjust='sidak')

methods = list(data.keys())
n_methods = len(methods)

alpha = 0.05

print("Comparing %s with %d runs each and alpha: %f" % (methods, len(data[methods[0]]), alpha))

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
            print("metod_a: %15s  method_b: %15s  p: %8.4f  mean_a: %6.2f  mean_b: %6.2f" % (
                method_a, method_b, p, mean_a, mean_b
            ))
