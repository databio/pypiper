""" Tests for effects of pipeline manager's halt() function. """

import os
from pypiper.flags import PAUSE_FLAG


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def test_halt_state(get_pipe_manager):
    """ Requesting a halt alters manager state. """
    pm = get_pipe_manager(name="test-pipe")
    assert pm.is_running
    assert pm._active
    pm.halt()
    assert pm.halted
    assert not pm.is_running
    assert not pm._active



def test_halt_file(get_pipe_manager):
    """ Requesting a halt produces a particular flag file. """
    pm = get_pipe_manager(name="TestPM")
    path_halt_file = pm.flag_file_path(PAUSE_FLAG)
    assert not os.path.isfile(path_halt_file)
    pm.halt()
    assert os.path.isfile(path_halt_file)
