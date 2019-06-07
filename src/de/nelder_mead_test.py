import numpy as np
import pytest
import sys
import os

from rosetta_pack import RosettaPack
from protein_data import ProteinData

from local_search.nelder_mead import nelder_mead

pname = '1crn'
rp = RosettaPack(name=pname)
eps = 1e-8


class TestProtein:
    def test_improves_once(self):
        pd = ProteinData(rp)

        score_before = pd(pd.angles)
        _, solution = nelder_mead(pd)
        score_after = pd(solution)

        assert score_after < score_before

    def test_improves_twice(self):
        pd = ProteinData(rp)

        score_before = pd(pd.angles)
        _, solution = nelder_mead(pd)
        score_after = pd(solution)
        _, solution = nelder_mead(pd)
        score_after_2 = pd(solution)

        assert score_after < score_before
        assert score_after_2 < score_after
