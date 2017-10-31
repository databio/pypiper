""" Tests for construction of checkpoint filepath """

import glob
import os
from pypiper import PipelineManager
from pypiper.stage import CHECKPOINT_EXTENSION, Stage
from tests.helpers import named_param


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



class DummyPM(PipelineManager):
    """ Simple override of true PipelineManager, for __init__ simplicity """
    def __init__(self, name, outfolder):
        self.name = name
        self.outfolder = outfolder



class CheckpointFilepathTests:
    """ Tests for determination of checkpoint filepath. """


    @named_param(argnames=["name1", "name2"],
                 argvalues=[("chipseq", "ATACseq"), ("rnaKallisto", "wgbs")])
    @named_param(argnames="spec_type",
                 argvalues=["stage_name", "stage", "function"])
    def test_distinguishes_pipelines_within_outfolder(
            self, name1, name2, spec_type, tmpdir):
        """ Checkpoint files within sample folder include pipeline name. """

        def trim_reads():
            pass

        def stage_spec():
            if spec_type == "function":
                return trim_reads
            elif spec_type not in ["stage", "stage_name"]:
                raise ValueError("Unrecognized stage specification type: {}".
                                 format(spec_type))
            else:
                s = Stage(trim_reads)
                return s.name if spec_type == "stage_name" else s

        outfolder = tmpdir.strpath

        # At start, we should have no checkpoints.
        checkpoint_pattern = os.path.join(outfolder, "*" + CHECKPOINT_EXTENSION)
        assert [] == glob.glob(checkpoint_pattern)

        plm1 = DummyPM(name1, outfolder)
        plm2 = DummyPM(name2, outfolder)

        checkpoint_name = "trim_reads"
        plm1.checkpoint(stage_spec())

        # Find the checkpoints; there should only be one.
        checkpoints = glob.glob(checkpoint_pattern)
        assert 1 == len(checkpoints)
        # Check that we have the expected checkpoint.
        exp_chkpt_fpath = os.path.join(outfolder, "{}_{}".format(
                name1, checkpoint_name + CHECKPOINT_EXTENSION))
        assert exp_chkpt_fpath == checkpoints[0]

        # Create a second checkpoint with the same stage, but
        plm2.checkpoint(stage_spec())
        checkpoints = glob.glob(checkpoint_pattern)
        assert 2 == len(checkpoints)
        exp_chkpt_fpath_2 = os.path.join(outfolder, "{}_{}".format(
                name2, checkpoint_name + CHECKPOINT_EXTENSION))
        assert {exp_chkpt_fpath, exp_chkpt_fpath_2} == set(checkpoints)