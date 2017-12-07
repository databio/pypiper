""" Tests for the timestamp functionality of a PipelineManager. """

import argparse
import mock
import os
import sys

import pytest

from pypiper import PipelineManager
from pypiper.utils import checkpoint_filepath, CHECKPOINT_SPECIFICATIONS


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
        """ If in particular state, managed pipeline halts on timestamp(). """
        pm = get_pipe_manager(name="TestPM")
        assert pm.is_running
        pm.halt_on_next = True
        pm.timestamp("testing")
        assert not pm.is_running
        assert pm.halted


    def test_correctly_sets_halt_on_next(self, get_pipe_manager):
        """ Of critical importance to timestamp's checkpointing functionality
        is its ability to alter the manager's state such that it triggers a
        halt on the subsequent timestamp() call. This allows timestamp() to
        be used in a prospective fashion while still preserving the ability to
        specify an inclusive rather than exclusive stopping point/stage.
        That is, it allows a pipeline to adopt the convention of calling
        timestamp() before beginning a conceptual block of processing logic,
        yet still (behave as though) stopping just after completion of
        execution of a defined stopping point. Essentially, the timestamp()
        calls can be prospective yet mixed with a retrospective halt point. """
        pm = get_pipe_manager(name="TestPM")
        pm.stop_after = "step2"
        assert pm.is_running
        assert not pm.halt_on_next
        pm.timestamp(checkpoint="step1")
        assert pm.is_running
        assert not pm.halt_on_next
        pm.timestamp(checkpoint="step2")
        assert pm.is_running
        assert pm.halt_on_next




class TimestampCheckpointTests:
    """ Tests for the checkpoint component of a timestamp() call. """
    pass



class TimestampStatusTypeTests:
    """ Tests for the type of status that a timestamp() call represents. """
    pass
