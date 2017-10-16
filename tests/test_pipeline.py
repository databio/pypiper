""" Tests for the Pipeline data type """

import os
import pytest
from pypiper import Pipeline


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


FILE1 = "file1.txt"
FILE2 = "file2.txt"


@pytest.fixture
def pl_mgr():
    pass



class DummyPipeline(Pipeline):
    """ Basic pipeline implementation for tests """

    def __init__(self, manager, output_folder):
        super(DummyPipeline, self).__init__("test-pipe", manager=manager)
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
    pass
