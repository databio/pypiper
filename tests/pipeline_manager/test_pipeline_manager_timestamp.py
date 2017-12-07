""" Tests for the timestamp functionality of a PipelineManager. """

import argparse
import mock
import os
import sys

import pytest

from pypiper import PipelineManager
from pypiper.utils import CHECKPOINT_SPECIFICATIONS


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def test_timestamp_requires_no_arguments(get_pipe_manager):
    """ A call to timestamp() requires no arguments. """
    pm = get_pipe_manager(name="TestPM")
    pm.timestamp()



def test_timestamp_message(get_pipe_manager):
    """ Tests for the message component of a timestamp() call. """
    pm = get_pipe_manager(name="TestPM")
    logfile = pm.pipeline_log_file
    with open(logfile, 'r') as f:
        orig_lines = f.readlines()
    pm.timestamp()



class TimestampHaltingTests:
    """ Tests for a manager's ability to halt a pipeline. """

    def test_halts_if_halt_on_next(self, get_pipe_manager):
        pm = get_pipe_manager(name="TestPM")
        assert pm.is_running
        pm.timestamp("testing")
        pm.halt_on_next = True




class TimestampCheckpointTests:
    """ Tests for the checkpoint component of a timestamp() call. """
    pass



class TimestampStatusTypeTests:
    """ Tests for the type of status that a timestamp() call represents. """
    pass
