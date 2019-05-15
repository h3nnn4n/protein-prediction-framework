import numpy as np
import pytest
import sys
import os

from rosetta_pack import RosettaPack
from protein_data import ProteinData

from local_search.hooke_jeeves import hooke, hooke_

pname = '1crn'
rp = RosettaPack(name=pname)
eps = 1e-8


class TestRosenbrock:
    def test_rosenbrock_improoves(self):
        start_point = np.array([-1.2, 1.0])
        nvars = 2
        initial_value = rosenbrock(start_point, nvars)

        itermax = 5000
        rho = 0.5
        eps = 1.0E-04

        _, solution = hooke_(nvars, start_point, rho, eps, itermax, rosenbrock)

        final_value = rosenbrock(solution, nvars)

        assert final_value < initial_value

    def test_rosenbrock_reaches_zero(self):
        start_point = np.array([-1.2, 1.0])
        nvars = 2

        itermax = 5000
        rho = 0.5
        eps = 1.0E-04

        _, solution = hooke_(nvars, start_point, rho, eps, itermax, rosenbrock)

        final_value = rosenbrock(solution, nvars)

        assert final_value < 1e-4


class TestProtein:
    def test_improves_once(self):
        pd = ProteinData(rp)

        score_before = pd(pd.angles)
        _, solution = hooke(pd)
        score_after = pd(solution)

        assert score_after < score_before

    def test_improves_twice(self):
        pd = ProteinData(rp)

        score_before = pd(pd.angles)
        _, solution = hooke(pd)
        score_after = pd(solution)
        _, solution = hooke(pd)
        score_after_2 = pd(solution)

        assert score_after < score_before
        assert score_after_2 < score_after


def rosenbrock(x, n):
    value = 100.0 * (x[1] - x[0] * x[0]) ** 2 \
        + (1.0 - x[0]) ** 2

    return value
