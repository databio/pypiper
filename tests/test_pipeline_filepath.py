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

    # Select the pipeline name.
    if "pipe_name" in request.fixturenames:
        pipe_name = request.getfixturevalue("pipe_name")
    else:
        pipe_name = "test-pipe"

    # Set output folder and name attributes for mocked PipelineManager.
    mock_mgr = mock.Mock(outfolder=tmpdir.strpath)
    type(mock_mgr).name = pipe_name  # Circumvent 'name' keyword on Mock.
    return mock_mgr


def test_requires_filename_or_suffix(pl_mgr):
    """Either filename or suffix is required to build a path."""
    with pytest.raises(TypeError):
        pipeline_filepath(pl_mgr)


@pytest.mark.parametrize(argnames="pipe_name", argvalues=PIPELINE_NAMES)
@pytest.mark.parametrize(argnames="suffix", argvalues=SUFFICES)
@pytest.mark.parametrize(
    argnames="test_type", argvalues=["has_pipe_name", "has_suffix", "full_path"]
)
def test_uses_pipeline_name_if_no_filename(
    pipe_name, suffix, test_type, pl_mgr, tmpdir
):
    """Pipeline name is proxy for filename if just suffix is given."""

    observed = pipeline_filepath(pl_mgr, suffix=suffix)

    # Allow test type to determine assertion.
    if test_type == "has_pipe_name":
        assert pipe_name in observed
    elif test_type == "has_suffix":
        assert observed.endswith(suffix)
    elif test_type == "full_path":
        try:
            expected = os.path.join(tmpdir.strpath, pipe_name + suffix)
            assert expected == observed
        except AssertionError:
            print("OUTFOLDER: {}".format(pl_mgr.outfolder))
            raise
    else:
        raise ValueError("Unrecognized test type: '{}'".format(test_type))


@pytest.mark.parametrize(
    argnames="filename", argvalues=["testfile" + suffix for suffix in SUFFICES]
)
@pytest.mark.parametrize(argnames="test_type", argvalues=["filename", "filepath"])
def test_direct_filename(tmpdir, filename, pl_mgr, test_type):
    """When given, filename is used instead of pipeline name."""
    fullpath = pipeline_filepath(pl_mgr, filename=filename)
    if test_type == "filename":
        _, observed = os.path.split(fullpath)
        assert filename == observed
    elif test_type == "filepath":
        expected = os.path.join(tmpdir.strpath, filename)
        assert expected == fullpath
    else:
        raise ValueError("Unrecognized test type: '{}'".format(test_type))


@pytest.mark.parametrize(argnames="filename", argvalues=["output", "testfile"])
@pytest.mark.parametrize(argnames="suffix", argvalues=SUFFICES)
def test_suffix_is_appended_to_filename_if_both_are_provided(pl_mgr, filename, suffix):
    """Suffix is appended to filename if both are provided."""
    expected = filename + suffix
    fullpath = pipeline_filepath(pl_mgr, filename=filename, suffix=suffix)
    _, observed = os.path.split(fullpath)
    assert expected == observed
