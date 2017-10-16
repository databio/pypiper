""" Tests for the Pipeline data type """

import os
import pytest
from pypiper import Pipeline, PipelineManager
from pypiper.manager import COMPLETE_FLAG


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
    pm = PipelineManager(name=TEST_PIPE_NAME, outfolder=tmpdir.strpath)
    def _ensure_stopped():
        pm.stop_pipeline()
    request.addfinalizer(_ensure_stopped)
    return pm



@pytest.fixture
def dummy_pipe(pl_mgr):
    """ Provide a basic Pipeline instance for a test case. """
    return DummyPipeline(pl_mgr)



class DummyPipeline(Pipeline):
    """ Basic pipeline implementation for tests """

    def __init__(self, manager):
        super(DummyPipeline, self).__init__(TEST_PIPE_NAME, manager=manager)

    def stages(self):
        return [self._write_file1, self._write_file2]

    def _write_file1(self):
        self._write(*FILE_TEXT_PAIRS[0])

    def _write_file2(self):
        self._write(*FILE_TEXT_PAIRS[1])

    def _write_file3(self):
        self._write(*FILE_TEXT_PAIRS[2])

    def _write(self, filename, content):
        path = os.path.join(self.manager.outfolder, filename)
        with open(path, 'w') as f:
            f.write(content)



class MostBasicPipelineTests:
    """ Test pipeline defined with notion of 'absolute minimum' config. """

    def test_runs_through_full(self, dummy_pipe, test_type, tmpdir):
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
            pass
        elif test_type == "pipe_flag":
            exp_flag = os.path.join(tmpdir.strpath, COMPLETE_FLAG)
            assert os.path
        else:
            raise ValueError("Unknown test type: {}".format(test_type))

    def test_skip_completed(self, dummy_pipe, test_type, tmpdir):
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
