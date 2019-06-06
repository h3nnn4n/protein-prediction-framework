# import os
# import re
# import sys

# import numpy as np
# import pandas as pd
# import seaborn as sns

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from utils import *


sns.set(style="whitegrid");


set_test('tests_qualify_with_ffi9_vs_old')

proteins = list_proteins()

experiments = ['sade_mc',
               'sade_remc',
               'sade_mc_ffi9_02',
               'sade_remc_ffi9_02']

for p in proteins:
    set_protein(p)
    
    for e in experiments:
        set_experiment(e)

        set_env()

        plot_all(metrics['best_rmsd'])
        
        if p in ['1crn', '1ail']:
            plt.ylim(0, 15)
        else:
            plt.ylim(0, 10)
        save_fig('best_rmsd')

        plot_all_best_and_mean_score_and_avg_rmsd()
