""" Tests for a pipeline's ability to checkpoint its stages. """

import os
import time

from pypiper.utils import checkpoint_filepath
from tests.helpers import fetch_checkpoint_files, named_param

from .conftest import get_pipeline

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


def test_pipeline_checkpoint_respect_sensitivity_checkpoint_perspective(
    pl_name, tmpdir
):
    """Pipeline can skip past its stage(s) for which checkpoint exists."""

    # Create the pipeline.
    pipeline = get_pipeline(pl_name, tmpdir.strpath)

    # Negative control to start test, that we have no checkpoint files.
    assert [] == fetch_checkpoint_files(pipeline.manager)

    # Generate some checkpoints.
    pipeline.run()

    # Verify that we created each of the checkpoints.
    expected = [
        checkpoint_filepath(f.__name__, pipeline.manager) for f in pipeline.functions
    ]
    observed = fetch_checkpoint_files(pipeline.manager)
    assert set(expected) == set(observed)

    # Collect checkpoint file timestamps for comparison after second run.
    timestamps = {f: os.path.getmtime(f) for f in observed}

    # Remove the checkpoint for the final stage.
    last_aligner_stage = pipeline.functions[-1]
    last_aligner_checkfile = checkpoint_filepath(last_aligner_stage, pipeline.manager)
    os.unlink(last_aligner_checkfile)

    # Verify removal of final stage checkpoint file.
    assert all([os.path.isfile(f) for f in expected[:-1]])
    assert not os.path.exists(last_aligner_checkfile)
    assert set(expected) != set(fetch_checkpoint_files(pipeline.manager))

    # Delay briefly so that we can more reliably compare checkpoint file
    # timestamps after a second pipeline run.
    time.sleep(0.05)

    # Repeat the pipeline's execution, but now with checkpoint file(s) for a
    # subset of its stages in place.
    pipeline.run()

    # Verify that we've restored the full collection of the pipeline's
    # checkpoint files to existence.
    observed = fetch_checkpoint_files(pipeline.manager)
    exp = set(expected)
    obs = set(observed)
    assert set(expected) == set(
        observed
    ), "Expected only:\n{}\nExpected and observed:\n{}\nObserved only:\n{}".format(
        exp - obs, exp & obs, obs - exp
    )

    # Verify the we didn't recreate the checkpoint file for each skipped stage.
    for f in expected[:-1]:
        expected_timestamp = timestamps[f]
        observed_timestamp = os.path.getmtime(f)
        assert expected_timestamp == observed_timestamp

    # Verify the we did in fact recreate the checkpoint file for the stage
    # that was rerun.
    assert (
        os.path.getmtime(last_aligner_checkfile) > timestamps[last_aligner_checkfile]
    ), "Recreated checkpoint file ('{}') should be newer than original".format(
        last_aligner_checkfile
    )


def test_pipeline_checkpoint_sensitivity_effect_perspective(pl_name, tmpdir):
    """The pipeline skips execution of stages with extant checkpoint."""

    # Create the pipeline, then check creation of output file.
    pipeline = get_pipeline(pl_name, tmpdir.strpath)
    output_file = os.path.join(pipeline.outfolder, pipeline.name_output_file)
    assert not os.path.exists(output_file)
    pipeline.run()
    assert os.path.isfile(output_file)

    # Validate pipeline effects (output file content).
    with open(output_file, "r") as f:
        lines = f.readlines()
    assert [s.name + os.linesep for s in pipeline.stages()] == lines

    # Verify presence of checkpoint files to support our expectation about
    # which stages should be skipped and which should be run during the second
    # time through the pipeline's execution.
    exp_cp_fpaths = set(
        checkpoint_filepath(s.name, pipeline.manager) for s in pipeline.stages()
    )
    assert exp_cp_fpaths == set(fetch_checkpoint_files(pipeline.manager))
    final_stage = pipeline.stages()[-1]
    final_stage_fpath = checkpoint_filepath(final_stage.name, pipeline.manager)
    os.unlink(final_stage_fpath)

    # Verify the effect of the second execution of the pipeline.
    pipeline.run()
    with open(output_file, "r") as f:
        lines = f.readlines()
    assert [final_stage.name + os.linesep] == lines


@named_param("overwrite", [False, True])
def test_pipeline_reruns_downstream_stages_according_to_parameterization(
    overwrite, pl_name, tmpdir
):
    """Pipeline overwrites downstream stages unless configured otherwise."""

    pl = get_pipeline(pl_name, tmpdir.strpath)

    # Create checkpoint file for each stage.
    stage_names = [s.name for s in pl.stages()]
    assert 1 < len(
        stage_names
    ), "Need pipeline with at least two stages to run this test."
    for s_name in stage_names:
        open(checkpoint_filepath(s_name, pl.manager), "w").close()

    # Remove the checkpoint file for the penultimate stage.
    penultimate_stage = stage_names[-2]
    os.unlink(checkpoint_filepath(penultimate_stage, pl.manager))

    # Configure the pipeline based on parameterization and run it starting
    # from the penultimate stage.
    pl.manager.overwrite_checkpoints = overwrite
    pl.run(start_point=penultimate_stage)

    # If we're overwriting downstream checkpoints, the last two stages are
    # run while otherwise only the penultimate stage is run.
    exp_stages = [stage_names[-2]]
    if overwrite:
        exp_stages.append(stage_names[-1])
    exp_lines = [func + os.linesep for func in stage_names[-2:]]
    outpath = os.path.join(pl.outfolder, pl.name_output_file)
    with open(outpath, "r") as f:
        obs_lines = f.readlines()
    assert exp_lines == obs_lines
