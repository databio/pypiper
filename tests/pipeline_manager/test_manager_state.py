""" Tests related to pipeline manager state. """

import os

import pytest

from pypiper.utils import checkpoint_filepath, pipeline_filepath
from tests.helpers import named_param

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


def test_starts_running(get_pipe_manager):
    """A PipelineManager begins running during its construction."""
    pm = get_pipe_manager(name="TestPM")
    assert pm._active


# Parameters governing execution:
# 1 -- checkpoint existence
# 3 -- halt state (.halted)


class ExecutionSkippingTests:
    """Tests for cases in which command execution should be skipped."""

    @named_param("start_point", ["align_reads", "make_call"])
    def test_skips_to_start(self, get_pipe_manager, start_point):
        """The pipeline manager can skip to a starting point."""

        # Initialize the manager.
        pm = get_pipe_manager(name="StartTestPM", start_point=start_point)

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
        cmds = [
            "fastqc",
            "touch {}".format(fastqc_rawfile),
            "touch {}".format(fastqc_zipfile),
        ]
        pm.run(cmds, target=fastqc_zipfile)
        assert not os.path.isfile(fastqc_zipfile)
        assert not os.path.isfile(fastqc_rawfile)

        # Make a all that should be the first executed, on the basis of
        # being associated with the designated.
        pm.timestamp(checkpoint=start_point)
        path_first_file = pipeline_filepath(pm, filename="outfile.bam")
        cmd = "touch {}".format(path_first_file)
        pm.run(cmd, target=path_first_file)
        assert os.path.isfile(path_first_file)

    @named_param("num_skips", argvalues=[1, 2, 3])
    def test_skips_execution_if_in_unstarted_state(self, get_pipe_manager, num_skips):
        """Pipeline manager skips command execution if not in active state."""

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

    @named_param("num_skips", argvalues=[1, 2, 3])
    def test_respects_checkpoints(self, get_pipe_manager, num_skips):
        """Manager can skip pipeline to where it's not yet checkpointed."""

        pm = get_pipe_manager(name="respect-checkpoints")

        # Control for possibility that skips are due to being in inactive mode.
        assert pm._active

        stages = ["merge", "qc", "filter", "align", "call"]

        # Create checkpoints.
        for s in stages[:num_skips]:
            pm.timestamp(checkpoint=s)

        # Go through the stages and see that we're skipping checkpoints
        # that exist, then proceeding to execute each subsequent stage.
        for i, s in enumerate(stages):
            outfile = pipeline_filepath(pm, s + ".txt")
            cmd = "touch {}".format(outfile)
            pm.timestamp(checkpoint=s)
            pm.run(cmd, target=outfile)

            if i < num_skips:
                # We should not have created the output file.
                try:
                    assert not os.path.isfile(outfile)
                except AssertionError:
                    print("Have run {} stage(s) of {} skip(s)".format(i + 1, num_skips))
                    print("Current manager checkpoint: {}".format(pm.curr_checkpoint))
                    raise
            else:
                # We should have created the output file.
                try:
                    assert os.path.isfile(outfile)
                except AssertionError:
                    print("Have run {} stage(s) of {} skip(s)".format(i + 1, num_skips))
                    print("Current manager checkpoint: {}".format(pm.curr_checkpoint))
                    print("Active? {}".format(pm._active))
                    raise

    @named_param("halt_index", [1, 2, 3])
    def test_respects_halt(self, get_pipe_manager, halt_index):
        """The pipeline manager skips execution if it's in halted state."""
        pm = get_pipe_manager(name="respects-halt")
        targets = ["file{}.txt".format(i) for i in range(1, 5)]
        for i, t in enumerate(targets):
            if i == halt_index:
                pm.halt(raise_error=False)
            target = pipeline_filepath(pm, filename=t)
            cmd = "touch {}".format(target)
            pm.run(cmd, target=target)
        for i, t in enumerate(targets):
            target = pipeline_filepath(pm, filename=t)
            if i < halt_index:
                assert os.path.isfile(target)
            else:
                assert not os.path.isfile(target)
