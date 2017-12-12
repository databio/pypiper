""" Tests for case in which multiple pipelines process a single sample. """

import pytest
from pypiper.utils import checkpoint_filepath
from tests.helpers import fetch_checkpoint_files
from .conftest import get_peak_caller, get_read_aligner


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def test_checkpoints_are_pipeline_unique(tmpdir):
    """ Names of checkpoint files depend on both stage and pipeline. """

    # Create two different pipelines.
    align_reads = get_read_aligner(tmpdir.strpath)
    call_peaks = get_peak_caller(tmpdir.strpath)

    # Get the stage names associated with each pipeline.
    alignment_stage_names = set(map(lambda s: s.name, align_reads.stages()))
    peak_call_stage_names = set(map(lambda s: s.name, call_peaks.stages()))

    # Check that we have one specific stage name shared between the pipelines.
    assert {"align_reads"} == alignment_stage_names & peak_call_stage_names
    assert align_reads.outfolder == call_peaks.outfolder

    # We begin with no checkpoint files.
    assert [] == list(fetch_checkpoint_files(align_reads.manager))
    assert [] == list(fetch_checkpoint_files(call_peaks.manager))

    # Run each pipeline.
    align_reads.run()
    call_peaks.run()

    # We expect a different checkpoint file for each stage of each pipeline.
    align_reads_expected = {checkpoint_filepath(s.name, align_reads)
                            for s in align_reads.stages()}
    call_peaks_expected = {checkpoint_filepath(s.name, call_peaks)
                           for s in call_peaks.stages()}

    # Pipeline names are unique here, and each checkpoint name includes
    # pipeline name for disambiguation, so even a pair of pipelines with a
    # nonempty stage name intersection has an empty checkpoint filenames
    # intersection, so long as the pipeline names are unique.
    assert set() == (align_reads_expected & call_peaks_expected)

    # When not setting start/stop parameters and beginning with no checkpoint
    # files in place, each pipeline generates its full set of checkpoint files.
    expected_checkpoints = align_reads_expected | call_peaks_expected
    observed_checkpoints = set(fetch_checkpoint_files(align_reads)) | \
                           set(fetch_checkpoint_files(call_peaks))

    # Verify satisfaction of expectation.
    try:
        assert expected_checkpoints == observed_checkpoints
    except AssertionError:
        only_exp = expected_checkpoints - observed_checkpoints
        exp_and_obs = expected_checkpoints & observed_checkpoints
        only_obs = observed_checkpoints - expected_checkpoints
        print("Only in expected:\n{}".format("\n".join(only_exp)))
        print("Expected and observed:\n{}".format("\n".join(exp_and_obs)))
        print("Only in observed:\n{}".format("\n".join(only_obs)))
        raise



@pytest.mark.skip("not implemented")
def test_pipeline_checkpoint_respect_specificity():
    pass
