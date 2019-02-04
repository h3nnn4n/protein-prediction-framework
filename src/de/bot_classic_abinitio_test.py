import pytest
import mock
from bot_classic_abinitio import main


def test_it_raises_when_no_arg():
    with pytest.raises(NotImplementedError):
        with mock.patch('sys.argv', ['main.py']):
            main()


def test_it_runs():
    setup_test_todo()
    with mock.patch('sys.argv', ['main.py', 'bot_rosetta_test_dummy_todo', 1]):
        with mock.patch('bot_classic_abinitio.ClassicAbinitio', get_de_boot_mock()):
            main()


# Utils

def get_de_boot_mock():
    boot_mock = mock.MagicMock(returns=None)
    return boot_mock


# Lazy poor mans fixture. FIXIT, someday
def setup_test_todo():
    with open('bot_rosetta_test_dummy_todo', 'wt') as f:
        f.write('1crn 1\n')
        f.write('1ail 1\n')
