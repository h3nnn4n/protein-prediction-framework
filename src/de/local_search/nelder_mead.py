from scipy.optimize import minimize


def nelder_mead(protein, eps=1e-04, max_evals=5000):
    disp = False

    x0 = protein.angles
    result = minimize(protein, x0, method='Nelder-Mead', tol=eps, options={'maxiter': max_evals, 'disp': disp})

    print(result['nfev'])

    return result['nfev'], protein.angles
