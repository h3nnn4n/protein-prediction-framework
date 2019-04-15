import os
import re
import sys

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from utils import *


base_path = os.getcwd()

sns.set(style="whitegrid")

set_experiment_path('tests_renan_vs_rosetta_12prots')
set_protein('1crn')
set_experiment('1crn_sade_remc')


_, mean_score = average_all(metrics['mean_score'])
_, best_score = average_all(metrics['best_score'])
evals, avg_rmsd = average_all(metrics['avg_rmsd'])
evals, moment_of_inertia = average_all(metrics['moment_of_inertia'])

plt.figure(figsize=(16.5, 10))

sns.lineplot(evals, mean_score, color='g')
sns.lineplot(evals, best_score, color='b')

ax2 = plt.twinx()

sns.lineplot(evals, avg_rmsd, color='r', ax=ax2)
# sns.lineplot(evals, best_score, label='best_score')

# plt.ylim(0, 125)
plt.xlim(0, 500000)


os.chdir(base_path)
plt.savefig('out.png')
