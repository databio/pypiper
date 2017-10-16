""" Tests for the Pipeline data type """

import os
import pytest
from pypiper import Pipeline, PipelineManager


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


FILE1 = "file1.txt"
FILE2 = "file2.txt"
TEST_PIPE_NAME = "test-pipe"



@pytest.fixture
def pl_mgr(request):
    """ Provide a PipelineManager and ensure that it's stopped. """
    pm = PipelineManager(TEST_PIPE_NAME)
    def _ensure_stopped():
        pass
    request.addfinalizer(_ensure_stopped)



class PipelineManagerTestContext(object):
    """ Ensure a test case's PipelineManager shuts down. """

    def __init__(self):
        pass

    def __enter__(self):
        pass

    def __exit__(self, *args, **kwargs):
        pass


class DummyPipeline(Pipeline):
    """ Basic pipeline implementation for tests """

    def __init__(self, manager, output_folder):
        super(DummyPipeline, self).__init__(TEST_PIPE_NAME, manager=manager)
        self.outpath = output_folder

    def stages(self):
        return [self._write_file1, self._write_file2]

    def _write_file1(self):
        self._write(FILE1, "hello there")

    def _write_file2(self):
        self._write(FILE2, "hello2")

    def _write(self, filename, content):
        path = os.path.join(self.outpath, filename)
        with open(path, 'w') as f:
            f.write(content)



class MostBasicPipelineTests:
    """ Test pipeline defined with notion of 'absolute minimum' config. """

    def test_runs_through_full(self):
        pass

    def test_skip_completed(self):
        pass

    def test_start_point(self):
        pass

    def test_start_before_completed_checkpoint(self):
        pass

    def test_same_start_stop(self):
        pass

    def test_stop_before_start(self):
        pass

    def test_can_skip_downstream_completed(self):
        pass

    def test_can_rerun_downstream_completed(self):
        pass

    def test_stop_at(self):
        pass

    def test_stop_after(self):
        pass
