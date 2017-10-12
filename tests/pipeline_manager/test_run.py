""" Tests for PipelineManager's run() method """

import pytest
from pypiper import PipelineManager
from pypiper.stage import CHECKPOINT_EXTENSION


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



@pytest.fixture
def pl_mgr(request, tmpdir):
    if "name" in request.fixturenames:
        pipe_name = request.getfixturevalue("name")
    else:
        pipe_name = "test-pipe"
    if "outfolder" in request.fixturenames:
        outfolder = request.getfixturevalue("outfolder")
    else:
        outfolder = tmpdir.strpath
    # Set 'multi' to prevent interference with stdout/err.
    pm = PipelineManager(pipe_name, outfolder, multi=True)
    return pm



class RunCheckpointTests:
    """ Tests for PipelineManager.run() with respect to checkpoint """


    @pytest.mark.parametrize(
        argnames="checkpoint", argvalues=["phase1", "stage2"])
    def test_checkpoint_name_is_satisfied(
            self, pl_mgr, checkpoint, tmpdir, capsys):
        """ The PipelineManager should respect checkpoint file. """



    def test_checkpoint_name_not_satisfied(self):
        """ Checkpoint file can be inferred from name.  """
        pass


    def test_checkpoint_file_is_satisfied(self):
        """ Execution stops if checkpoint file exists. """
        pass


    def test_checkpoint_file_not_satisfied(self):
        """ Execution proceeds if the checkpoint file doesn't exist. """
        pass


    def test_checkpoint_file_precludes_name(self):
        """ The PipelineManager should continue if there's no checkpoint file. """
        pass


    def test_no_checkpoint_specified(self):
        """ Execution proceeds if command to run is not checkpointed. """
        pass
