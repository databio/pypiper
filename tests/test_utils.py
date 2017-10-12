""" Tests for utility functions """

import os
import mock
import pytest
from pypiper.utils import pipeline_filepath


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


PIPELINE_NAMES = ["chiapet", "chipseq", "atacseq", "kallisto", "wgbs"]
SUFFICES = [".txt", "_results.csv", ".stats.tsv", "-data.json"]



@pytest.fixture
def pl_mgr(request, tmpdir):
    """
    Provide test case with a mocked PipelineManager instance.

    :param pytest.fixtures.FixtureRequet request: test case requesting the
        setup fixture / parameterization
    :param py.path.local.LocalPath tmpdir: Test case temporary path object.
    :return mock.MagicMock: Mocked PipelineManager, sufficient for test.
    """

    # Select the pipline name.
    if "pipe_name" in request.fixturenames:
        pipe_name = request.getfixturevalue("pipe_name")
    else:
        pipe_name = "test-pipe"

    # Set output folder and name attributes for mocked PipelineManager.
    mock_mgr = mock.Mock(outfolder=tmpdir.strpath)
    type(mock_mgr).name = pipe_name    # Circumvent 'name' keyword on Mock.
    return mock_mgr



class PipelineFilepathTests:
    """ Test creation of filepath within managed pipeline's output folder. """


    def test_requires_filename_or_suffix(self, pl_mgr):
        """ Either filename or suffix is required to build a path. """
        with pytest.raises(TypeError):
            pipeline_filepath(pl_mgr)


    @pytest.mark.parametrize(argnames="pipe_name", argvalues=PIPELINE_NAMES)
    @pytest.mark.parametrize(argnames="suffix", argvalues=SUFFICES)
    def test_uses_pipeline_name_if_no_filename(
            self, pipe_name, suffix, pl_mgr, tmpdir):
        """ Pipeline name is proxy for filename if just suffix is given. """
        expected = os.path.join(tmpdir.strpath, pipe_name + suffix)
        observed = pipeline_filepath(pl_mgr, suffix=suffix)
        try:
            assert expected == observed
        except AssertionError:
            print("OUTFOLDER: {}".format(pl_mgr.outfolder))
            raise


    def test_suffix_if_no_filename(self, pl_mgr):
        """ Suffix is appended to pipeline name if filename isn't given. """
        pass


    def test_direct_filename(self, pl_mgr):
        """ When given, filename is used instead of pipeline name. """
        pass


    def test_suffix_is_appended_to_filename_if_both_are_provided(self, pl_mgr):
        """ Suffix is appended to filename if both are provided. """
        pass