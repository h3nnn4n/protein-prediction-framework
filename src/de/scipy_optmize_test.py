import pickle
from scipy.optimize import minimize, dual_annealing
from protein_data import ProteinData
from rosetta_pack import RosettaPack


attributes_to_save = ['fun', 'nfev', 'nit']


def run_minimize(n_runs=40):
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
                results[method] = {}

            for attr in attributes_to_save:
                if not hasattr(result, attr):
                    continue

                if attr not in results[method].keys():
                    results[method][attr] = []

                value = getattr(result, attr)

                if method == 'Powell' and attr == 'fun':
                    value = float(value)

                results[method][attr].append(value)

            print("%20s %8d %8.3f" % (method, result['nfev'], result['fun']))

    return results


def run_others(n_runs=40):
    rp = RosettaPack(name='1zdd')
    pd = ProteinData(rp)

    results = {}

    run_dual_annealing(results=results, pd=pd, n_runs=n_runs)

    return results


def run_dual_annealing(results=None, pd=None, n_runs=40):
    method = 'Dual Annealing'

    maxiter = 20
    maxfun = 5000

    print('Running ' + method)
    for _ in range(n_runs):
        pd.reset()
        x0 = pd.angles
        bounds = [(-180, 180) for _ in range(len(x0))]
        result = dual_annealing(pd, bounds=bounds, x0=x0, maxiter=maxiter, maxfun=maxfun)

        if method not in results.keys():
            results[method] = {}

        for attr in attributes_to_save:
            if not hasattr(result, attr):
                continue

            if attr not in results[method].keys():
                results[method][attr] = []

            value = getattr(result, attr)

            if method == 'Powell' and attr == 'fun':
                value = float(value)

            results[method][attr].append(value)

        print("%20s %8d %8.3f" % (method, result['nfev'], result['fun']))

    return results


data = run_others(n_runs=2)

with open('scipy_min_dual_annealing.pickle', 'wb') as f:
    pickle.dump(data, f)

data = run_minimize(n_runs=10)

with open('scipy_min_test.pickle', 'wb') as f:
    pickle.dump(data, f)
