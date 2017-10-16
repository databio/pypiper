""" Tests for the Pipeline data type """

from functools import partial
import os
import pytest
from pypiper import Pipeline, PipelineManager
from pypiper.manager import COMPLETE_FLAG
from pypiper.pipeline import checkpoint_filepath
from pypiper.stage import Stage
from pypiper.utils import flag_name


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


FILE1_NAME = "file1.txt"
FILE2_NAME = "file2.txt"
FILE3_NAME = "file3.txt"
FILENAMES = [FILE1_NAME, FILE2_NAME, FILE3_NAME]

FILE1_TEXT = "hello there"
FILE2_TEXT = "hello2"
FILE3_TEXT = "third"
CONTENTS = [FILE1_TEXT, FILE2_TEXT, FILE3_TEXT]

FILE_TEXT_PAIRS = list(zip(FILENAMES, CONTENTS))
TEST_PIPE_NAME = "test-pipe"


def pytest_generate_tests(metafunc):
    if "test_type" in metafunc.fixturenames and \
            metafunc.cls == MostBasicPipelineTests:
        metafunc.parametrize(argnames="test_type", 
                             argvalues=["effects", "checkpoints", "pipe_flag"])


@pytest.fixture
def pl_mgr(request, tmpdir):
    """ Provide a PipelineManager and ensure that it's stopped. """
    pm = PipelineManager(
            name=TEST_PIPE_NAME, outfolder=tmpdir.strpath, multi=True)
    def _ensure_stopped():
        pm.stop_pipeline()
    request.addfinalizer(_ensure_stopped)
    return pm


@pytest.fixture
def dummy_pipe(pl_mgr):
    """ Provide a basic Pipeline instance for a test case. """
    return DummyPipeline(pl_mgr)


def _write_file1(folder):
    _write(*FILE_TEXT_PAIRS[0], folder=folder)


def _write_file2(folder):
    _write(*FILE_TEXT_PAIRS[1], folder=folder)


def _write_file3(folder):
    _write(*FILE_TEXT_PAIRS[2], folder=folder)


def _write(filename, content, folder=None):
    path = os.path.join(folder, filename)
    with open(path, 'w') as f:
        f.write(content)


class DummyPipeline(Pipeline):
    """ Basic pipeline implementation for tests """

    def __init__(self, manager):
        super(DummyPipeline, self).__init__(TEST_PIPE_NAME, manager=manager)

    @property
    def stages(self):
        # File content writers parameterized with output folder.
        fixed_folder_funcs = []
        for f in [_write_file1, _write_file2, _write_file3]:
            f_fixed = partial(f, folder=self.outfolder)
            f_fixed.__name__ = f.__name__
            fixed_folder_funcs.append(f_fixed)
        return [Stage(f) for f in fixed_folder_funcs]



class MostBasicPipelineTests:
    """ Test pipeline defined with notion of 'absolute minimum' config. """


    def test_runs_through_full(self, dummy_pipe, test_type, tmpdir):
        dummy_pipe.run(start=None, stop_at=None, stop_after=None)
        if test_type == "effects":
            fpath_text_pairs = [(os.path.join(tmpdir.strpath, fname), content)
                                for fname, content in FILE_TEXT_PAIRS]
            for fpath, content in fpath_text_pairs:
                assert os.path.isfile(fpath)
                exp_content = content.split(os.linesep)
                with open(fpath, 'r') as f:
                    obs_content = [l.rstrip(os.linesep) for l in f.readlines()]
                assert exp_content == obs_content
        elif test_type == "checkpoints":
            for stage in [_write_file1, _write_file2, _write_file3]:
                chkpt_fpath = checkpoint_filepath(stage, dummy_pipe)
                try:
                    assert os.path.isfile(chkpt_fpath)
                except AssertionError:
                    print("Stage '{}' file doesn't exist: '{}'".format(
                        stage.__name__, chkpt_fpath))
                    raise
        elif test_type == "pipe_flag":
            flags = os.listdir(dummy_pipe.outfolder)
            assert 1 == len(flags)
            exp_flag = os.path.join(tmpdir.strpath, flag_name(COMPLETE_FLAG))
            assert os.path.isfile(exp_flag)
        else:
            raise ValueError("Unknown test type: {}".format(test_type))


    def test_skip_completed(self, dummy_pipe, test_type, tmpdir):
        """ Pre-completed stage(s) are skipped. """
        pass

    def test_start_point(self, dummy_pipe, test_type, tmpdir):
        pass

    def test_start_before_completed_checkpoint(self, dummy_pipe, test_type, tmpdir):
        pass

    def test_same_start_stop(self, dummy_pipe, test_type, tmpdir):
        pass

    def test_stop_before_start(self, dummy_pipe, test_type, tmpdir):
        pass

    def test_can_skip_downstream_completed(self, dummy_pipe, test_type, tmpdir):
        pass

    def test_can_rerun_downstream_completed(self, dummy_pipe, test_type, tmpdir):
        pass

    def test_stop_at(self, dummy_pipe, test_type, tmpdir):
        pass

    def test_stop_after(self, dummy_pipe, test_type, tmpdir):
        pass
