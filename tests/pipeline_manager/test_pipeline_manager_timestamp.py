""" Tests for the timestamp functionality of a PipelineManager. """

import os
import sys
import pytest
from pypiper.utils import checkpoint_filepath
from tests.helpers import fetch_checkpoint_files, named_param


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


STATE_TEST = "state"
FILES_TEST = "files"



def pytest_generate_tests(metafunc):
    """ Dynamic test case generation for this module. """
    if "retrospective" in metafunc.fixturenames:
        metafunc.parametrize("retrospective", [False, True])
    if "test_type" in metafunc.fixturenames:
        metafunc.parametrize("test_type", [FILES_TEST, STATE_TEST])



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
    # the luxury of a logfile for the pipeline manager.
    out, err = capsys.readouterr()
    sys.stdout.write(out)
    sys.stderr.write(err)

    # The stdout capture with capsys comes through as a single unicode block.
    assert message_content in str(out), \
            "Missing timestamp message ('{}') in message(s)".\
            format(message_content)



class TimestampHaltingTests:
    """ Tests for a manager's ability to halt a pipeline. """


    # Note that the tests here are not truly logically independent from the
    # functionality of the manager's halt() method. The assertions made here
    # assume a particular effects of a call to the halt() method, and that
    # that method is behaving in accordance with expectation. For true
    # independence, we could mock halt() and make assertions about calls on
    # the mock, but here that seems to inject a level of complexity for which
    # the cost exceeds the benefit of the logical independence that it confers.



    def test_halts_if_hitting_exclusive_halt_point(self, get_pipe_manager):
        """ Halt point may be specified prospectively. """

        # Create manager, set halt point, and check that it's running.
        halt_name = "phase3"
        pm = get_pipe_manager(name="TestPM")
        pm.stop_before = halt_name
        assert pm.is_running

        # Make non-halting checkpointed timestamp calls.
        pm.timestamp(checkpoint="phase1")
        assert pm.is_running
        assert not pm.halted
        pm.timestamp(checkpoint="phase2")
        assert pm.is_running
        assert not pm.halted

        # Make the halt-inducing checkpointed timestamp call.
        pm.timestamp(checkpoint=halt_name)

        # Verify that we've halted.
        try:
            assert pm.halted
        except AssertionError:
            print("STATUS: {}".format(pm.status))
            raise
        assert not pm.is_running


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

        # Establish manager and perform initial control assertions.
        pm = get_pipe_manager(name="TestPM")
        pm.stop_after = "step2"
        assert pm.is_running
        assert not pm.halt_on_next

        # Make non-halt-status-altering checkpointed timestamp call and
        # verify that we're still running and that we're not scheduled to halt.
        pm.timestamp(checkpoint="step1")
        assert pm.is_running
        assert not pm.halt_on_next

        # Make halt-status-altering checkpointed timestamp call and verify
        # that we're still running and that we've now been scheduled to halt.
        pm.timestamp(checkpoint="step2")
        assert pm.is_running
        assert pm.halt_on_next



class TimestampStatusTypeTests:
    """ Tests for the type of status that a timestamp() call represents. """


    def test_initial_timestamp_checkpoint_file(
            self, get_pipe_manager, retrospective):
        """ Initial checkpointed timestamp writes checkpoint file if and only
        if it's a retrospective timestamp. """
        pm = get_pipe_manager(name="init-timestamp-file")
        stage_name = "align_reads"
        pm.timestamp(checkpoint=stage_name, finished=retrospective)
        check_fpath = checkpoint_filepath(stage_name, pm)
        if retrospective:
            assert os.path.isfile(check_fpath)
        else:
            assert not os.path.isfile(check_fpath)


    @named_param("which_checkpoint_state",
                 ["curr_checkpoint", "prev_checkpoint"])
    def test_initial_timestamp_states(
            self, get_pipe_manager, retrospective, which_checkpoint_state):
        """ Which checkpoint state is updated by a checkpointed timestamp
        call depends upon the perspective of the call. """

        # Create the manager and make the timestamp call.
        pm = get_pipe_manager(name="InitialTimestampState")
        stage_name = "quality_control"
        pm.timestamp(checkpoint=stage_name, finished=retrospective)

        # Form expectations.
        if retrospective:
            prev_exp = stage_name
            curr_exp = None
        else:
            prev_exp = None
            curr_exp = stage_name

        # Make the assertion on the specified checkpoint state.
        if which_checkpoint_state == "curr_checkpoint":
            assert curr_exp == getattr(pm, "curr_checkpoint")
        else:
            assert prev_exp == getattr(pm, "prev_checkpoint")


    def test_two_prospective_checkpointed_timestamps(
            self, test_type, stage_pair, pm):
        """ Prospective timestamp generates file for previous checkpoint. """

        stage1, stage2 = stage_pair
        pm.timestamp(checkpoint=stage1, finished=False)
        pm.timestamp(checkpoint=stage2, finished=False)

        if test_type == FILES_TEST:
            checkpoint_files = fetch_checkpoint_files(pm)
            expected = [checkpoint_filepath(stage1, pm)]
            assert set(expected) == set(checkpoint_files)
        else:
            assert stage1 == pm.prev_checkpoint
            assert stage2 == pm.curr_checkpoint


    def test_two_retrospective_checkpointed_timestamps(
            self, test_type, stage_pair, pm):
        """ Retrospective timestamp generates file for current checkpoint. """

        stage1, stage2 = stage_pair
        pm.timestamp(checkpoint=stage1, finished=True)
        pm.timestamp(checkpoint=stage2, finished=True)

        if test_type == FILES_TEST:
            checkpoint_files = fetch_checkpoint_files(pm)
            expected = [checkpoint_filepath(s, pm) for s in [stage1, stage2]]
            assert set(expected) == set(checkpoint_files)
        else:
            assert stage2 == pm.prev_checkpoint
            assert pm.curr_checkpoint is None


    def test_prospective_then_retrospective_checkpointed_timestamps(
            self, test_type, stage_pair, pm):
        """ If a prospective checkpointed timestamp is followed by a
        retrospective one, there's only a file for the retrospective one. """

        stage1, stage2 = stage_pair
        pm.timestamp(checkpoint=stage1, finished=False)
        assert stage1 == pm.curr_checkpoint
        pm.timestamp(checkpoint=stage2, finished=True)

        if test_type == FILES_TEST:
            checkpoint_files = fetch_checkpoint_files(pm)
            expected = [checkpoint_filepath(stage2, pm)]
            assert set(expected) == set(checkpoint_files)
        else:
            # Current checkpoint will be reset by second (retrospective)
            # timestamp call.
            assert stage2 == pm.prev_checkpoint
            assert pm.curr_checkpoint is None


    @pytest.mark.skip("not implemented")
    def test_retrospective_the_prospective_checkpointed_timestamps(
            self, test_type, stage_pair, pm):
        stage1, stage2 = stage_pair
        pm.timestamp(checkpoint=stage1, finished=True)
        pm.timestamp(checkpoint=stage2, finished=False)
        if test_type == FILES_TEST:
            assert os.path.isfile(checkpoint_filepath(stage1, pm))
            assert os.path.isfile(checkpoint_filepath(stage2, pm))
        else:
            assert stage1 == pm.prev_checkpoint
            assert stage2 == pm.curr_checkpoint


    @pytest.fixture
    def stage_pair(self):
        return "merge_input", "quality_control"


    @pytest.fixture
    def pm(self, get_pipe_manager):
        return get_pipe_manager(name="checkpointed-timestamp-pair")
