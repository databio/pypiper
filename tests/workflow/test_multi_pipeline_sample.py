""" Tests for case in which multiple pipelines process a single sample. """

import os
import time
import pytest
from pypiper import Stage
from pypiper.utils import checkpoint_filepath
from tests.helpers import fetch_checkpoint_files, SafeTestPipeline


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


READ_ALIGNER_FILENAME = "aligner.lst"
PEAK_CALLER_FILENAME = "caller.lst"



def pytest_generate_tests(metafunc):
    """ Dynamic test case parameterization for this module. """
    if "pl_name" in metafunc.fixturenames:
        metafunc.parametrize(
            "pl_name", [read_aligner.__name__, call_peaks.__name__])



# Dummy functions used as elements of pipeline stages() collections.
def merge_input():
    pass

def qc():
    pass

def align_reads():
    pass

def call_peaks():
    pass



class FunctionNameWriterPipeline(SafeTestPipeline):
    """ Basic pipeline that writes to file the names of its functions. """
    
    
    def __init__(self, name, outfolder, filename, functions):
        """
        Name and outfolder go to generic pipeline ctor; filename and functions
        collection are used specifically by instances of this class.

        :param str name: Name for this pipeline.
        :param str outfolder: Path to pipeline's output folder.
        :param str filename: Name for file in which to write function names.
        :param Sequence[callable] functions: Functions on which this pipeline
            is to operate (i.e., the functions for which name should be
            written to output file).
        """
        # Set instance-specific variables.
        self.name_output_file = filename
        self.functions = functions
        # Get the stages() benefit of superclass extension.
        super(FunctionNameWriterPipeline, self).__init__(
                name=name, outfolder=outfolder)

    def write_name(self, func):
        """
        Write the name of a function to this pipeline's output file.

        :param callable func: Name of function to write to the output file.
        """
        outpath = os.path.join(self.outfolder, self.name_output_file)
        with open(outpath, 'a') as f:
            f.write(func.__name__ + os.linesep)

    def run(self):
        """ Start with clean output file, then use superclass method. """
        # Ensure that we start with a clean file since the nature of the
        # operations performed (sequential file writes) creates desire to
        # open output file in append mode rather than write mode.
        output_file = os.path.join(self.outfolder, self.name_output_file)
        if os.path.exists(output_file):
            os.unlink(output_file)
        super(FunctionNameWriterPipeline, self).run()


    def stages(self):
        """ Sequence of operations to perform. """
        return [Stage(self.write_name, (f, ), name=f.__name__)
                for f in self.functions]


# Functions and fixtures

def _get_read_aligner(outfolder):
    """ Create a dummy 'read aligner' pipeline. """
    return FunctionNameWriterPipeline(
            "read-aligner", outfolder,
            READ_ALIGNER_FILENAME, [merge_input, qc, align_reads])


def _get_peak_caller(outfolder):
    """ Create a dummy 'peak caller' pipeline. """
    return FunctionNameWriterPipeline(
        "peak-caller", outfolder,
        PEAK_CALLER_FILENAME, [align_reads, call_peaks])

def _get_pipeline(name, outfolder):
    if name == read_aligner.__name__:
        return _get_read_aligner(outfolder)
    elif name == call_peaks.__name__:
        return _get_peak_caller(outfolder)
    else:
        raise ValueError("Unknown pipeline request: '{}'".format(name))


@pytest.fixture
def read_aligner(tmpdir):
    """ Provide test case with a read aligner pipeline instance. """
    return _get_read_aligner(outfolder=tmpdir.strpath)


@pytest.fixture
def peak_caller(tmpdir):
    """ Provide test case with a 'PeakCaller' pipeline instance. """
    return _get_peak_caller(outfolder=tmpdir.strpath)



def test_pipeline_checkpoint_respect_sensitivity_checkpoint_perspective(
        pl_name, tmpdir):
    """ Pipeline can skip past its stage(s) for which checkpoint exists. """

    # Create the pipeline.
    pipeline = _get_pipeline(pl_name, tmpdir.strpath)

    # Negative control to start test, that we have no checkpoint files.
    assert [] == fetch_checkpoint_files(pipeline.manager)

    # Generate some checkpoints.
    pipeline.run()

    # Verify that we created each of the checkpoints.
    expected = [checkpoint_filepath(f.__name__, pipeline.manager)
                for f in pipeline.functions]
    observed = fetch_checkpoint_files(pipeline.manager)
    assert set(expected) == set(observed)

    # Collect checkpoint file timestamps for comparison after second run.
    timestamps = {f: os.path.getmtime(f) for f in observed}

    # Remove the checkpoint for the final stage.
    last_aligner_stage = pipeline.functions[-1]
    last_aligner_checkfile = checkpoint_filepath(
            last_aligner_stage, pipeline.manager)
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
    assert set(expected) == set(observed), \
            "Expected only:\n{}\nExpected and observed:\n{}\nObserved only:\n{}".format(
                    exp - obs, exp & obs, obs - exp)

    # Verify the we didn't recreate the checkpoint file for each skipped stage.
    for f in expected[:-1]:
        expected_timestamp = timestamps[f]
        observed_timestamp = os.path.getmtime(f)
        assert expected_timestamp == observed_timestamp

    # Verify the we did in fact recreate the checkpoint file for the stage
    # that was rerun.
    assert os.path.getmtime(last_aligner_checkfile) > \
           timestamps[last_aligner_checkfile], \
            "Recreated checkpoint file ('{}') should be newer than original".\
           format(last_aligner_checkfile)



def test_pipeline_checkpoint_sensitivity_effect_perspective(pl_name, tmpdir):
    """ The pipeline skips execution of stages with extant checkpoint. """

    # Create the pipeline, then check creation of output file.
    pipeline = _get_pipeline(pl_name, tmpdir.strpath)
    output_file = os.path.join(pipeline.outfolder, pipeline.name_output_file)
    assert not os.path.exists(output_file)
    pipeline.run()
    assert os.path.isfile(output_file)

    # Validate pipeline effects (output file content).
    with open(output_file, 'r') as f:
        lines = f.readlines()
    assert [s.name + os.linesep for s in pipeline.stages()] == lines

    # Verify presence of checkpoint files to support our expectation about
    # which stages should be skipped and which should be run during the second
    # time through the pipeline's execution.
    exp_cp_fpaths = set(checkpoint_filepath(s.name, pipeline.manager)
                        for s in pipeline.stages())
    assert exp_cp_fpaths == set(fetch_checkpoint_files(pipeline.manager))
    final_stage = pipeline.stages()[-1]
    final_stage_fpath = checkpoint_filepath(final_stage.name, pipeline.manager)
    os.unlink(final_stage_fpath)

    # Verify the effect of the second execution of the pipeline.
    pipeline.run()
    with open(output_file, 'r') as f:
        lines = f.readlines()
    assert [final_stage.name + os.linesep] == lines



@pytest.mark.skip("not implemented")
def test_pipeline_can_overwrite_checkpoints():
    pass



@pytest.mark.skip("not implemented")
def test_different_pipeline_checkpoints_are_unique_for_the_sample(tmpdir):
    pass



@pytest.mark.skip("not implemented")
def test_pipeline_checkpoint_respect_specificity():
    pass
