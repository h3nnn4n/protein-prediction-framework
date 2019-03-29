import pytest
import mock

from forced_insertion import ForcedInsertion
from rosetta_pack import RosettaPack
from protein_data import ProteinData

pname = '1crn'
rp = RosettaPack(name=pname)
eps = 1e-8


class TestInit:
    def test_no_params(self):
        ForcedInsertion()

    def test_with_params(self):
        mock_de = mock.MagicMock
        mock_de.rosetta_pack = rp

        ForcedInsertion(
            de=mock_de,
            mode='9mer',
            chance=0.1,
            enable=True
        )

    def test_raises_for_unknown_mode(self):
        mock_de = mock.MagicMock
        mock_de.rosetta_pack = rp

        with pytest.raises(ValueError):
            ForcedInsertion(
                de=mock_de,
                mode='potato',
                chance=0.1,
                enable=True
            )


class TestCanRun:
    def test_should_run(self):
        mock_de = mock.MagicMock
        mock_de.rosetta_pack = rp

        fi = ForcedInsertion(
            de=mock_de,
            mode='9mer',
            chance=0.1,
            enable=True
        )

        assert fi.can_run()

    def test_zero_chance(self):
        mock_de = mock.MagicMock
        mock_de.rosetta_pack = rp

        fi = ForcedInsertion(
            de=mock_de,
            mode='9mer',
            chance=0.0,
            enable=True
        )

        assert not fi.can_run()

    def test_no_de(self):
        fi = ForcedInsertion(
            mode='9mer',
            chance=0.1,
            enable=True
        )

        assert not fi.can_run()

    def test_no_mode(self):
        mock_de = mock.MagicMock
        mock_de.rosetta_pack = rp

        fi = ForcedInsertion(
            de=mock_de,
            chance=0.1,
            enable=True
        )

        assert not fi.can_run()

    def test_enable_override(self):
        mock_de = mock.MagicMock
        mock_de.rosetta_pack = rp

        fi = ForcedInsertion(
            de=mock_de,
            mode='9mer',
            chance=0.1,
            enable=False
        )

        assert not fi.can_run()


class TestRun:
    def test_runs(self):
        mock_de = de_mock_builder()

        fi = ForcedInsertion(
            de=mock_de,
            mode='9mer',
            chance=0.1,
            enable=True
        )

        fi.run()

    def test_it_changes_the_confirmation_with_9mer(self):
        mock_de = de_mock_builder()

        fi = ForcedInsertion(
            de=mock_de,
            mode='9mer',
            chance=0.1,
            enable=True
        )

        fi.logger_file_handle = mock.MagicMock
        fi.logger_file_handle.write = mock.MagicMock(return_value=True)
        fi.logger_file_handle.flush = mock.MagicMock(return_value=True)

        for _ in range(10):
            old_angles = rp.get_backbone_angles_from_pose(mock_de.pop[0].pose)

            with mock.patch("forced_insertion.uniform", mock.MagicMock(return_value=0.01)):
                fi.run()

            new_angles = rp.get_backbone_angles_from_pose(mock_de.pop[0].pose)

            different_angles = 0

            for a, b in zip(new_angles, old_angles):
                if abs(a - b) > eps:
                    different_angles += 1

            assert different_angles == 9 * 3


def pop_builder(pop_size=1):
    return [ProteinData(rp) for _ in range(pop_size)]


def de_mock_builder():
    mock_de = mock.MagicMock
    mock_de.rosetta_pack = rp
    mock_de.pop = pop_builder()
    mock_de.spent_evals = 0
    mock_de.it = 0

    return mock_de
