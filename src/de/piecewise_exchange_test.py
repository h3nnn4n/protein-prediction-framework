import pytest
import mock

from protein_data import ProteinData
from rosetta_pack import RosettaPack
from piecewise_exchange import PiecewiseExchange

pname = '1crn'
rp = RosettaPack(name=pname)
ss_pred = rp.ss_pred


class TestSplit:
    def test_split3(self):
        mock_de = self.de_mock_builder()
        mock_de.rosetta_pack.ss_pred = 'CCCHHHCCCHHHCCC'

        pe = PiecewiseExchange(de=mock_de)
        splits = pe.split_regions()

        assert splits == [3, 6, 9, 12]

    def test_split4(self):
        mock_de = self.de_mock_builder()
        mock_de.rosetta_pack.ss_pred = 'CCCCHHHHCCCCHHHHCCCC'

        pe = PiecewiseExchange(de=mock_de)
        splits = pe.split_regions()

        assert splits == [4, 8, 12, 16]

    def de_mock_builder(self):
        mock_de = mock.MagicMock
        mock_de.rosetta_pack = rp

        return mock_de


class TestSwapRegion:
    def test_swap_first_region(self):
        mock_de = self.de_mock_builder()

        pe = PiecewiseExchange(de=mock_de)
        splits = pe.split_regions()

        for i in range(splits[0]):
            assert mock_de.pop[1].angles[i] == 10
            assert mock_de.pop[5].angles[i] == 50

        pe.swap_region([0] + splits[0:1], 1, 5)

        for i in range(splits[0]):
            assert mock_de.pop[1].angles[i] == 50
            assert mock_de.pop[5].angles[i] == 10

    def test_swap_second_region(self):
        mock_de = self.de_mock_builder()

        pe = PiecewiseExchange(de=mock_de)
        splits = pe.split_regions()

        for i in range(splits[0], splits[1]):
            assert mock_de.pop[1].angles[i] == 10
            assert mock_de.pop[5].angles[i] == 50

        pe.swap_region(splits[0:2], 1, 5)

        for i in range(splits[0], splits[1]):
            assert mock_de.pop[1].angles[i] == 50
            assert mock_de.pop[5].angles[i] == 10

    def de_mock_builder(self):
        mock_de = mock.MagicMock
        mock_de.rosetta_pack = rp
        mock_de.pop = self.build_pop()
        mock_de.rosetta_pack.ss_pred = 'CCCHHHCCC'

        return mock_de

    def build_pop(self):
        pop_size = 10

        pop = [ProteinData(rp) for _ in range(pop_size)]

        for i in range(len(rp.target)):
            for j in range(pop_size):
                pop[j].angles[i] = j * 10

        return pop


class TestScramblePop:
    def test_scramble_pop(self):
        mock_de = self.de_mock_builder()

        pe = PiecewiseExchange(de=mock_de)
        _ = [p.eval() for p in mock_de.pop]

        scores_before = [p.score for p in mock_de.pop]

        pe.scramble_pop()

        scores_after = [p.score for p in mock_de.pop]

        for a, b in zip(scores_before, scores_after):
            assert a != b

    def de_mock_builder(self):
        mock_de = mock.MagicMock
        mock_de.rosetta_pack = rp
        mock_de.pop = self.build_pop()
        mock_de.pop_size = len(mock_de.pop)
        mock_de.rosetta_pack.ss_pred = 'CCCHHHCCC'

        return mock_de

    def build_pop(self):
        pop_size = 10

        pop = [ProteinData(rp) for _ in range(pop_size)]

        for i in range(len(rp.target)):
            for j in range(pop_size):
                pop[j].angles[i] = j * 10

        return pop


class TestRandomPiecewiseSearch:
    def test_score_changed(self):
        mock_de = self.de_mock_builder()
        pe = PiecewiseExchange(de=mock_de)

        trial_score_before = mock_de.trial.score

        pe.random_piecewise_search()

        trial_score_after = mock_de.trial.score

        assert trial_score_before != trial_score_after

    def test_changed_all_angles(self):
        mock_de = self.de_mock_builder()
        pe = PiecewiseExchange(de=mock_de)

        pe.random_piecewise_search()

        for k, angle in enumerate(mock_de.trial.angles):
            assert angle != -50, "failed at angle %4d %8.3f" % (k, angle)

    def de_mock_builder(self):
        mock_de = mock.MagicMock
        mock_de.rosetta_pack = rp
        mock_de.rosetta_pack.ss_pred = ss_pred
        mock_de.pop = self.build_pop()
        mock_de.pop_size = len(mock_de.pop)
        mock_de.trial = self.build_trial()

        return mock_de

    def build_pop(self):
        pop_size = 10

        pop = [ProteinData(rp) for _ in range(pop_size)]

        for i in range(len(rp.target)):
            for j in range(pop_size):
                pop[j].angles[i] = j * 10

        return pop

    def build_trial(self):
        trial = ProteinData(rp)

        for i in range(len(rp.target)):
            trial.angles[i] = -50

        trial.eval()

        return trial
