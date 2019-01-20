import pickle

methods = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'L-BFGS-B', 'TNC']

with open('scipy_min_test.pickle', 'rb') as f:
    data = pickle.load(f)

print(data.keys())
