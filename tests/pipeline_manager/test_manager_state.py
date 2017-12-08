""" Tests related to pipeline manager state. """

import os
import pytest
from pypiper.utils import pipeline_filepath
from tests.helpers import named_param


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def test_starts_running(get_pipe_manager):
    """ A PipelineManager begins running during its construction. """
    pm = get_pipe_manager(name="TestPM")
    assert pm._active
    assert pm.is_running



@named_param("num_skips", argvalues=[1, 2, 3])
def test_skips_execution_if_in_unstarted_state(get_pipe_manager, num_skips):
    """ Pipeline manager skips command execution if not in active state. """

    pm = get_pipe_manager(name="skip-execs")
    pm._active = False

    testfile = pipeline_filepath(pm, filename="output.txt")
    assert not os.path.isfile(testfile)

    cmd = "touch {}".format(testfile)
    num_calls = 0

    # Remain inactive for a parameterized number of call-skipping iterations,
    # then adopt active mode.
    while True:
        pm.run(cmd, target=testfile)
        num_calls += 1
        if num_calls == num_skips:
            pm._active = True
        elif num_calls > num_skips:
            break
        # If we're still looping, we've not yet made a call in active mode.
        assert not os.path.isfile(testfile)

    # We break the loop once we've made a call in active state.
    assert os.path.isfile(testfile)


# Parameters governing execution:
# 1 -- checkpoint existence
# 2 -- start state (._active)
# 3 -- halt state (.halted)


class ExecutionSkippingTests:
    """ Tests for cases in which command execution should be skipped. """


    def test_skips_to_start(self, get_pipe_manager):
        """ The pipeline manager can skip to a starting point. """

        # Initialize the manager.
        pm = get_pipe_manager(name="StartTestPM", start_point="align_reads")

        # Make a call that should be skipped on the basis of not yet
        # reaching the start point.
        pm.timestamp(checkpoint="merge_reads")
        path_merge_file = pipeline_filepath(pm, filename="merge.txt")
        assert not os.path.isfile(path_merge_file)
        cmd = "touch {}".format(path_merge_file)
        pm.run(cmd, target=path_merge_file)
        assert not os.path.isfile(path_merge_file)

        # Make a call that should also be skipped on the basis of not yet
        # reaching the designated starting/activation point.
        pm.timestamp(checkpoint="fastqc")
        fastqc_folder = os.path.join(pm.outfolder, "fastqc")
        os.makedirs(fastqc_folder)
        fastqc_zipfile = os.path.join(fastqc_folder, "qc.zip")
        fastqc_rawfile = os.path.join(fastqc_folder, "qc.txt")
        cmds = ["fastqc", "touch {}".format(fastqc_rawfile),
                "touch {}".format(fastqc_zipfile)]
        pm.run(cmds, target=fastqc_zipfile)
        assert not os.path.isfile(fastqc_zipfile)
        assert not os.path.isfile(fastqc_rawfile)

        # Make a all that should be the first executed, on the basis of
        # being associated with the designated.
        pm.timestamp(checkpoint="align_reads")
        path_align_reads_file = pipeline_filepath(pm, filename="alignment.bam")
        cmd = "touch {}".format(path_align_reads_file)
        pm.run(cmd, target=path_align_reads_file)
        assert os.path.isfile(path_align_reads_file)


    @pytest.mark.skip("Not implemented")
    def test_respects_checkpoints(self, get_pipe_manager):
        """ Manager can skip pipeline to where it's not yet checkpointed. """
        pass


    @pytest.mark.skip("Not implemented")
    def test_respects_halt(self, get_pipe_manager):
        pass


    @pytest.mark.skip("Not implemented")
    def test_activation(self):
        pass
