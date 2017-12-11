""" Test configuration for Pipeline tests. """

import os
import pytest
from pypiper import Stage
from tests.helpers import SafeTestPipeline


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



READ_ALIGNER_FILENAME = "aligner.lst"
PEAK_CALLER_FILENAME = "caller.lst"



def pytest_generate_tests(metafunc):
    """ Dynamic test case parameterization. """
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
        return [Stage(self.write_name, (f,), name=f.__name__)
                for f in self.functions]



# Functions and fixtures

def get_read_aligner(outfolder):
    """ Create a dummy 'read aligner' pipeline. """
    return FunctionNameWriterPipeline(
        "read-aligner", outfolder,
        READ_ALIGNER_FILENAME, [merge_input, qc, align_reads])



def get_peak_caller(outfolder):
    """ Create a dummy 'peak caller' pipeline. """
    return FunctionNameWriterPipeline(
        "peak-caller", outfolder,
        PEAK_CALLER_FILENAME, [align_reads, call_peaks])



def get_pipeline(name, outfolder):
    """
    Build and return pipeline instance associated with given name.

    :param str name: Name of the pipeline to build.
    :param str outfolder: Path to output folder for use by pipeline instance.
    :return SafeTestPipeline: A test-session-safe instance of a Pipeline.
    """
    if name == read_aligner.__name__:
        return get_read_aligner(outfolder)
    elif name == call_peaks.__name__:
        return get_peak_caller(outfolder)
    else:
        raise ValueError("Unknown pipeline request: '{}'".format(name))



@pytest.fixture
def read_aligner(tmpdir):
    """ Provide test case with a read aligner pipeline instance. """
    return get_read_aligner(outfolder=tmpdir.strpath)



@pytest.fixture
def peak_caller(tmpdir):
    """ Provide test case with a 'PeakCaller' pipeline instance. """
    return get_peak_caller(outfolder=tmpdir.strpath)
