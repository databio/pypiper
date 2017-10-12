""" Tests for PipelineManager's run() method """

import os
import random
import string
import sys

import pytest

from pypiper import PipelineManager
from pypiper.stage import checkpoint_filename, checkpoint_filepath


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



class PipelineManagerTester(PipelineManager):
    """ When a test completes, ensure that a pipeline is stopped normally. """

    def run(self, *args, **kwargs):
        super(PipelineManagerTester, self).run(*args, **kwargs)
        self.stop_pipeline()



@pytest.fixture
def pl_mgr(request, tmpdir):
    """
    Provide a test case with a PipelineManager that cleanly stops.

    :param pytest.fixtures.FixtureRequest request: test case requesting this
        fixture parameterization
    :param py.path.local.LocalPath tmpdir: temporary test folder fixture
    :return PipelineManagerTester: thin wrapper around normal PipelineManager,
        that simply stops cleanly upon test cessation.
    """

    # Allow requesting test case to provide arguments for PipelineManager,
    # but also provide defaults.
    if "name" in request.fixturenames:
        pipe_name = request.getfixturevalue("name")
    else:
        pipe_name = "test-pipe"
    if "outfolder" in request.fixturenames:
        outfolder = request.getfixturevalue("outfolder")
    else:
        outfolder = tmpdir.strpath

    # Set 'multi' to prevent interference with stdout/err.
    pm = PipelineManagerTester(pipe_name, outfolder, multi=True)
    return pm



class RunCheckpointTests:
    """ Tests for PipelineManager.run() with respect to checkpoint """


    def _assert_checkpoint(self, captured, filepath=None):
        """ Check that an existing file message was captured. """
        message = "Checkpoint file exists"
        assert captured.startswith(message)
        if filepath is not None:
            assert filepath in captured


    @pytest.mark.parametrize(
        argnames="checkpoint", argvalues=["phase1", "stage2"])
    def test_checkpoint_name_is_satisfied(
            self, pl_mgr, checkpoint, tmpdir, capsys):
        """ The PipelineManager should respect checkpoint file. """

        # Create empty checkpoint file.
        open(tmpdir.join(checkpoint_filename(checkpoint)).strpath, 'w').close()

        # Create a dummy filepath that should not exist, but would if the
        # command did in fact run.
        path_dummy_file = tmpdir.join("dummy.txt").strpath
        assert not os.path.exists(path_dummy_file)

        # Make the call under test and ensure that we didn't create the
        # dummy file.
        pl_mgr.run("touch {}".format(path_dummy_file),
                   target=path_dummy_file, checkpoint=checkpoint)
        assert not os.path.exists(path_dummy_file)

        # Store captured output, then re-write for display in case of
        # test failure.
        out, err = capsys.readouterr()
        sys.stdout.write(out)
        sys.stdout.write(err)

        # Check that the expected message was generated.
        self._assert_checkpoint(captured=out)


    @pytest.mark.parametrize(
        argnames=["checkpoint", "previous"],
        argvalues=[("step2", ["part1"]),
                   ("phase_three", ["stage1", "part_zwei"])])
    def test_checkpoint_name_not_satisfied(
            self, checkpoint, previous, tmpdir, pl_mgr):
        """ Checkpoint file can be inferred from name.  """

        # Arbitrarily establish previous checkpoint files. This is irrelevant
        # for determination of whether a single run() call should execute,
        # but it implicitly tests specificity of checkpoint file name and
        # more accurately reflects a true use case.
        for prev in previous:
            name = checkpoint_filename(prev)
            path = tmpdir.join(name).strpath
            open(path, 'w').close()

        base_dummy_filename = "".join(
                [random.choice(string.ascii_lowercase) for _ in range(20)])
        path_dummy_file = tmpdir.join(base_dummy_filename + ".txt").strpath

        assert not os.path.exists(path_dummy_file)
        chkpt_file = checkpoint_filepath(checkpoint, pm=pl_mgr)
        pl_mgr.run("touch {}".format(path_dummy_file),
                   target=chkpt_file, checkpoint=checkpoint)

        assert os.path.exists(path_dummy_file)


    @pytest.mark.parametrize(
        argnames="chkpt_fname",
        argvalues=["start.completed", "middle.done", "final.finished"])
    @pytest.mark.parametrize(
        argnames="test_type", argvalues=["test_effect", "test_message"])
    def test_checkpoint_file_is_satisfied(
            self, chkpt_fname, pl_mgr, test_type, tmpdir, capsys):
        """ Execution stops if checkpoint file exists. """

        # Create command target and ensure it doesn't exist.
        dummy_target = tmpdir.join("dummy.txt").strpath
        assert not os.path.exists(dummy_target)

        # Create the empty checkpoint file.
        chkpt_fpath = tmpdir.join(chkpt_fname).strpath
        open(chkpt_fpath, 'w').close()

        # Make the call under test.
        cmd = "touch {}".format(dummy_target)
        pl_mgr.run(cmd, target=dummy_target, checkpoint_filename=chkpt_fname)

        # Conduct the test.
        if test_type == "test_effect":
            assert not os.path.exists(tmpdir.join(dummy_target).strpath)
        elif test_type == "test_message":
            out, err = capsys.readouterr()
            sys.stdout.write(out)
            sys.stderr.write(err)
            self._assert_checkpoint(captured=out, filepath=chkpt_fpath)
        else:
            raise ValueError("Unrecognized test type: '{}'".format(test_type))


    def test_checkpoint_file_not_satisfied(self):
        """ Execution proceeds if the checkpoint file doesn't exist. """
        pass


    def test_checkpoint_file_precludes_name(self):
        """ Direct checkpoint file takes precedence over name-derived ones. """
        pass


    def test_no_checkpoint_specified(self):
        """ Execution proceeds if command to run is not checkpointed. """
        pass
