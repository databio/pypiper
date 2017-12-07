""" Tests related to pipeline manager state. """

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def test_starts_running(get_pipe_manager):
    """ A PipelineManager begins running during its construction. """
    pm = get_pipe_manager(name="TestPM")
    assert pm.is_running
