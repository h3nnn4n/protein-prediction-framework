import pytest
import mock
from bot import main


def test_it_raises_when_no_arg():
    with pytest.raises(NotImplementedError):
        with mock.patch('sys.argv', ['main.py']):
            main()


def test_it_runs():
    setup_test_todo()
    with mock.patch('sys.argv', ['main.py', 'bot_test_dummy_todo', 1]):
        with mock.patch('bot.boot', get_de_boot_mock()):
            main()


# Utils

def get_de_boot_mock():
    boot_mock = mock.MagicMock(returns=None)
    return boot_mock


# Lazy poor mans fixture. FIXIT, someday
def setup_test_todo():
    with open('bot_test_dummy_todo', 'wt') as f:
        f.write('a.yaml\n')
        f.write('b.yaml\n')

    with open('a.yaml', 'wt') as f:
        f.write('hello world of tests\n')
