""" Tests for the Pipeline data type """

import pytest
from pypiper import Pipeline


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



class DummyPipeline(Pipeline):
    """ Basic pipeline implementation for tests """


    def stages(self):
        pass



class MostBasicPipelineTests:
    """ Test pipeline defined with notion of 'absolute minimum' config. """
    pass
