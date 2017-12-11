""" Tests for case in which multiple pipelines process a single sample. """

import os
import time
import pytest
from pypiper.utils import checkpoint_filepath
from tests.helpers import fetch_checkpoint_files
from .conftest import get_pipeline


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



@pytest.mark.skip("not implemented")
def test_different_pipeline_checkpoints_are_unique_for_the_sample(tmpdir):
    pass



@pytest.mark.skip("not implemented")
def test_pipeline_checkpoint_respect_specificity():
    pass
