import pytest
import mock

from crowding import Crowding
from rosetta_pack import RosettaPack
from protein_data import ProteinData

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


def pop_builder(pop_size=1):
    return [ProteinData(rp) for _ in range(pop_size)]


def de_mock_builder(pop_size=1):
    mock_de = mock.MagicMock
    mock_de.rosetta_pack = rp
    mock_de.pop = pop_builder(pop_size=pop_size)
    mock_de.pop_size = pop_size
    mock_de.spent_evals = 0
    mock_de.it = 0
    mock_de.trial = ProteinData(rp)

    return mock_de


def add_crowding_params(de, run_crowding=True, crowding_factor=5, crowding_mode='rmsd'):
    de.run_crowding = run_crowding
    de.crowding_factor = crowding_factor
    de.crowding_mode = crowding_mode
