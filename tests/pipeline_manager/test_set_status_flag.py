""" Tests for changes to pipepline manager's status flag. """

import pytest

from pypiper.flags import *
from pypiper.flags import __all__ as ALL_FLAGS
from tests.helpers import named_param


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



@named_param(argnames="status", argvalues=ALL_FLAGS)
def test_set_status_flag_is_idempotent(get_pipe_manager, status):
    """ Calls to manager's status flag setter are idempotent. """
    pm = get_pipe_manager(name="TestPM")
    pm._set_status_flag(status)
    assert status == pm.status
    pm._set_status_flag(status)
    assert status == pm.status



@pytest.mark.parametrize(
    argnames=["init_state", "new_state"],
    argvalues=[(WAIT_FLAG, RUN_FLAG), (WAIT_FLAG, COMPLETE_FLAG),
               (WAIT_FLAG, FAIL_FLAG), (RUN_FLAG, COMPLETE_FLAG),
               (RUN_FLAG, PAUSE_FLAG), (RUN_FLAG, FAIL_FLAG),
               (FAIL_FLAG, RUN_FLAG)])
def test_changes_status_state(get_pipe_manager, init_state, new_state):
    """ Manager setting status flag changes is internal status/state. """
    pm = get_pipe_manager(name="test-pipe")
    assert pm.status == RUN_FLAG
    pm._set_status_flag(init_state)
    assert init_state == pm.status
    pm._set_status_flag(new_state)
    assert new_state == pm.status
