import numpy as np
import pytest
import sys
import os

from random import random

from rosetta_pack import RosettaPack
from protein_data import ProteinData

from local_search.hooke_jeeves import hooke, hooke_

pname = '1crn'
rp = RosettaPack(name=pname)

pd = ProteinData(rp)

for i in range(10):
    rho = random()
    score_before = pd(pd.angles)
    evals, solution = hooke(pd, eps=1e-4, rho=rho)
    score_after = pd(solution)

    print('%3d %14.6f %14.6f %8d %8.8f' % (i, score_before, score_after, evals, rho))
