import pytest
import mock

from crowding import Crowding
from rosetta_pack import RosettaPack
from protein_data import ProteinData
from operators.operators import Operators

pname = '1zdd'
rp = RosettaPack(name=pname)


class TestCandidateSorting:
    def test_find_nearest_candidate_to_trial(self):
        de = de_mock_builder(pop_size=25)
        add_crowding_params(de, crowding_factor=5)
        trial = de.trial
        crowding = Crowding(de=de)

        candidates = crowding.get_random_individuals()
        nearest = crowding.find_nearest_candidate_to_trial(candidates)

        nearest_rmsd = rp.get_rmsd_from_pose(nearest.pose, trial.pose)

        for p in candidates:
            assert nearest_rmsd <= rp.get_rmsd_from_pose(p.pose, trial.pose)

    def test_sort_candidates_by_trial_rmsd(self):
        de = de_mock_builder(pop_size=25)
        add_crowding_params(de, crowding_factor=10)
        crowding = Crowding(de=de)
        trial = de.trial

        candidates = crowding.get_random_individuals()
        sorted_candidates = crowding.sort_candidates_by_trial_rmsd(candidates)

        for i in range(len(sorted_candidates) - 1):
            rmsd1 = rp.get_rmsd_from_pose(sorted_candidates[i].pose, trial.pose)
            rmsd2 = rp.get_rmsd_from_pose(sorted_candidates[i + 1].pose, trial.pose)
            assert rmsd1 <= rmsd2


class TestSwap:
    def test_didnt_duplicate(self):
        de = de_mock_builder(pop_size=6)
        add_crowding_params(de, crowding_factor=3)
        crowding = Crowding(de=de)
        trial = de.trial

        trial.stage2_mc()

        candidates = crowding.get_random_individuals()
        nearest = crowding.find_nearest_candidate_to_trial(candidates)

        crowding.swap_candidate_and_trial(nearest)

        trial.eval()
        nearest.eval()

        assert nearest.score != trial.score
        assert nearest.pose is not trial.pose

    def test_pop_consistency(self):
        de = de_mock_builder(pop_size=6)
        de.operators = Operators(de=de)

        op = de.operators.get_operator('monte_carlo_3s')

        for _ in range(10):
            for i in range(len(de.pop)):
                op(i)

        for i in range(len(de.pop)):
            for j in range(i + 1, len(de.pop)):
                assert i != j

                p1 = de.pop[i]
                p2 = de.pop[j]

                assert p1 is not p2
                assert p1 is not de.trial
                assert p2 is not de.trial

                assert p1.pose is not p2.pose
                assert p1.pose is not de.trial.pose
                assert p2.pose is not de.trial.pose


def pop_builder(pop_size=1):
    return [ProteinData(rp) for _ in range(pop_size)]


def de_mock_builder(pop_size=1):
    mock_de = mock.MagicMock
    mock_de.rosetta_pack = rp
    mock_de.pop = pop_builder(pop_size=pop_size)
    mock_de.pop_size = pop_size
    mock_de.spent_evals = 0
    mock_de.it = 0
    mock_de.sade_run = False
    mock_de.get_f_cr = lambda: (0.5, 0.5)
    add_crowding_params(mock_de, crowding_factor=3)
    mock_de.crowding = Crowding(de=mock_de)

    def selection_mock(candidate):
        mock_de.crowding.selection(candidate)

    mock_de.selection = selection_mock
    mock_de.trial = ProteinData(rp)

    return mock_de


def add_crowding_params(de, run_crowding=True, crowding_factor=5, crowding_mode='rmsd'):
    de.run_crowding = run_crowding
    de.crowding_factor = crowding_factor
    de.crowding_mode = crowding_mode
