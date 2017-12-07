""" Tests for the timestamp functionality of a PipelineManager. """

import argparse

import pytest

from pypiper import PipelineManager
from pypiper.utils import CHECKPOINT_SPECIFICATIONS


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



@pytest.fixture
def pl_mgr(request, tmpdir):
    """ Provide test case with a basic PipelineManager. """
    opts = {cp_spec: request.getfixturevalue(cp_spec) for cp_spec in
            CHECKPOINT_SPECIFICATIONS if cp_spec in request.fixturenames}
    args = argparse.Namespace(**opts)
    return PipelineManager("test-PM", outfolder=tmpdir.strpath, args=args)



def test_timestamp_requires_no_arguments():
    """ A call to timestamp() requires no arguments. """
    pass



class TimestampMessageTests:
    """ Tests for the message component of a timestamp() call. """
    pass



class TimestampCheckpointTests:
    """ Tests for the checkpoint component of a timestamp() call. """
    pass



class TimestampStatusTypeTests:
    """ Tests for the type of status that a timestamp() call represents. """
    pass
