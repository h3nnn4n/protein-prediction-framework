import pickle
from scipy.optimize import minimize
from protein_data import ProteinData
from rosetta_pack import RosettaPack


def run(n_runs=40):
    rp = RosettaPack(name='1zdd')
    pd = ProteinData(rp)

    tol = 1e-4

    methods = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'L-BFGS-B', 'TNC']

    maxiter = 10000
    disp = False

    results = {}

    for method in methods:
        print('Running ' + method)
        for _ in range(n_runs):
            pd.reset()
            x0 = pd.angles
            result = minimize(pd, x0, method=method, tol=tol, options={'maxiter': maxiter, 'disp': disp})

            if method not in results.keys():
                results[method] = []

            results[method].append(result.fun)
            print(result.fun)

    return results

data = run()

with open('scipy_min_test.pickle', 'wb') as f:
    pickle.dump(data, f)
