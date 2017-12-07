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



def test_timestamp_message(get_pipe_manager, capsys):
    """ Tests for the message component of a timestamp() call. """
    name = "TestPM"
    pm = get_pipe_manager(name=name)
    logfile = pm.pipeline_log_file
    assert not os.path.exists(logfile)
    message_content = "Just testing"
    message = "### {}".format(message_content)
    pm.timestamp(message)

    # Capture output but also write it so it's there.
    # Since we're in interactive mode for the testing session, we don't have
    # the luxury of a lofgile for the pipeline manager.
    out, err = capsys.readouterr()
    sys.stdout.write(out)
    sys.stderr.write(err)

    # The stdout capture with capsys comes through as a single unicode block.
    assert message_content in str(out), \
            "Missing timestamp message ('{}') in message(s)".\
            format(message_content)



class TimestampHaltingTests:
    """ Tests for a manager's ability to halt a pipeline. """

    def test_halts_if_halt_on_next(self, get_pipe_manager):
        pm = get_pipe_manager(name="TestPM")
        assert pm.is_running
        pm.halt_on_next = True
        pm.timestamp("testing")
        assert not pm.is_running
        assert pm.halted



class TimestampCheckpointTests:
    """ Tests for the checkpoint component of a timestamp() call. """
    pass



class TimestampStatusTypeTests:
    """ Tests for the type of status that a timestamp() call represents. """
    pass
