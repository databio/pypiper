""" Tests for effects of pipeline manager's halt() function. """

import os
import pytest
from pypiper.exceptions import PipelineHalt
from pypiper.flags import COMPLETE_FLAG, PAUSE_FLAG
from tests.helpers import named_param


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def test_halt_state(get_pipe_manager):
    """ Requesting a halt alters manager state. """
    pm = get_pipe_manager(name="test-pipe")
    assert pm._active
    pm.halt(raise_error=False)
    assert pm.halted
    assert not pm._active



def test_halt_file(get_pipe_manager):
    """ Requesting a halt produces a particular flag file. """
    pm = get_pipe_manager(name="TestPM")
    path_halt_file = pm._flag_file_path(PAUSE_FLAG)
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



@named_param("raise_error", [False, True])
@named_param("test_type", argvalues=["halt_flag", "complete_flag"])
def test_halt_status_supersedes_completed(
        get_pipe_manager, raise_error, test_type):
    """ Halting pipeline replaces completed flag with halt flag. """

    # Create manager and completion flag.
    pm = get_pipe_manager(name="halt-status-flag")
    pm._set_status_flag(COMPLETE_FLAG)
    path_complete_flag = pm._flag_file_path(COMPLETE_FLAG)
    assert os.path.isfile(path_complete_flag)

    # Perform the halt.
    try:
        pm.halt(raise_error=raise_error)
    except PipelineHalt:
        # We don't care about exceptionality here, just that the flag files
        # are adjusted regardless of the halt type.
        pass

    # Check either the presence of the halt flag or the absence of the
    # completion flag, depending on test parameterization.
    if test_type == "halt_flag":
        path_halt_flag = pm._flag_file_path(PAUSE_FLAG)
        assert os.path.isfile(path_halt_flag)
    elif test_type == "complete_flag":
        assert not os.path.isfile(path_complete_flag)
    else:
        raise ValueError("Unknown test type: '{}'".format(test_type))
