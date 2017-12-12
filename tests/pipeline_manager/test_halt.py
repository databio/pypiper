""" Tests for effects of pipeline manager's halt() function. """

import os
import pytest
from pypiper.exceptions import PipelineHalt
from pypiper.flags import PAUSE_FLAG
from tests.helpers import named_param


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def test_halt_state(get_pipe_manager):
    """ Requesting a halt alters manager state. """
    pm = get_pipe_manager(name="test-pipe")
    assert pm.is_running
    assert pm._active
    pm.halt(raise_error=False)
    assert pm.halted
    assert not pm.is_running
    assert not pm._active



def test_halt_file(get_pipe_manager):
    """ Requesting a halt produces a particular flag file. """
    pm = get_pipe_manager(name="TestPM")
    path_halt_file = pm.flag_file_path(PAUSE_FLAG)
    assert not os.path.isfile(path_halt_file)
    pm.halt(raise_error=False)
    assert os.path.isfile(path_halt_file)



@named_param("raise_error", [False, True, None])
def test_halt_exceptionality(get_pipe_manager, raise_error):
    """ Halting is conditionally exceptional """
    pm = get_pipe_manager(name="halt-error")
    if raise_error is None:
        # Default is exceptional.
        with pytest.raises(PipelineHalt):
            pm.halt()
    elif raise_error:
        with pytest.raises(PipelineHalt):
            pm.halt(raise_error=True)
    else:
        pm.halt(raise_error=False)

