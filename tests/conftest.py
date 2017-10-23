""" Fixtures and configuration visible to all tests """

from functools import partial
import os
import pytest
from pypiper import Pipeline, PipelineManager, Stage


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


# Use a weird suffix for glob specificity.
OUTPUT_SUFFIX = ".testout"

TEST_PIPE_NAME = "test-pipe"

FILE1_TEXT = "hello there"
FILE2_TEXT = "hello2"
FILE3_TEXT = "third"
CONTENTS = [FILE1_TEXT, FILE2_TEXT, FILE3_TEXT]

FILE1_NAME = "file1{}".format(OUTPUT_SUFFIX)
FILE2_NAME = "file2{}".format(OUTPUT_SUFFIX)
FILE3_NAME = "file3{}".format(OUTPUT_SUFFIX)
FILENAMES = [FILE1_NAME, FILE2_NAME, FILE3_NAME]

FILE_TEXT_PAIRS = list(zip(FILENAMES, CONTENTS))



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



def write_file1(folder):
    _write(*FILE_TEXT_PAIRS[0], folder=folder)



def write_file2(folder):
    _write(*FILE_TEXT_PAIRS[1], folder=folder)



def write_file3(folder):
    _write(*FILE_TEXT_PAIRS[2], folder=folder)



def _write(filename, content, folder=None):
    path = os.path.join(folder, filename)
    with open(path, 'w') as f:
        f.write(content)



class DummyPipeline(Pipeline):
    """ Basic pipeline implementation for tests """

    def __init__(self, manager):
        super(DummyPipeline, self).__init__(TEST_PIPE_NAME, manager=manager)

    def stages(self):
        """
        Establish the stages/phases for this test pipeline.

        :return list[pypiper.Stage]: Sequence of stages for this pipeline.
        """
        # File content writers parameterized with output folder.
        fixed_folder_funcs = []
        for f in [write_file1, write_file2, write_file3]:
            f_fixed = partial(f, folder=self.outfolder)
            f_fixed.__name__ = f.__name__
            fixed_folder_funcs.append(f_fixed)
        return [Stage(f) for f in fixed_folder_funcs]
