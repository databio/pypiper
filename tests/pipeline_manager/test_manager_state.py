""" Tests related to pipeline manager state. """

import os
import pytest
from pypiper.utils import pipeline_filepath


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def test_starts_running(get_pipe_manager):
    """ A PipelineManager begins running during its construction. """
    pm = get_pipe_manager(name="TestPM")
    assert pm.is_running


# Parameters governing execution:
# 1 -- checkpoint existence
# 2 -- start state (._has_started)
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
