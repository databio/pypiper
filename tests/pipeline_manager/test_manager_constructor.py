""" Test effects of construction of a pipeline manager. """

import argparse
import pytest
from tests.helpers import named_param


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def pytest_generate_tests(metafunc):
    """ Dynamic test case generation for this module's test cases. """
    if "spec_type" in metafunc.fixturenames:
        metafunc.parametrize(
                argnames="spec_type", argvalues=["cmdl", "kwargs"])



class ManagerConstructorCheckpointSpecificationTests:
    """ Tests for manager's constructor's ability to parse and set
    checkpoint specifications, which can determine aspects of control flow. """


    def test_no_checkpoint_specifications(self, get_pipe_manager):
        """ A manager may be constructed without any checkpoint provision. """
        get_pipe_manager(name="test-pipe")


    @named_param("start_point", ["filter_reads", "align_reads"])
    def test_just_start(self, get_pipe_manager, spec_type, start_point):
        """ Starting point may be set from command-line or ctor keyword. """
        spec_data = {"start_point": start_point}
        if spec_type == "cmdl":
            kwargs = {"args": argparse.Namespace(**spec_data)}
        else:
            kwargs = spec_data
        pm = get_pipe_manager(name="start-test", **kwargs)
        assert start_point == pm.start_point


    @named_param("stop_type", ["stop_before", "stop_after"])
    @named_param("stop_point", ["align_reads", "call_peaks"])
    def test_just_stop(self, get_pipe_manager,
                       spec_type, stop_type, stop_point):
        """ Particular stopping type is set correctly. """
        spec_data = {stop_type: stop_point}
        if spec_type == "cmdl":
            kwargs = {"args": argparse.Namespace(**spec_data)}
        else:
            kwargs = spec_data
        pm = get_pipe_manager(name="stop-test", **kwargs)
        assert stop_point == getattr(pm, stop_type)


    @pytest.mark.skip("Not implemented")
    def test_start_and_stop(self, get_pipe_manager, spec_type):
        """ Specifying both start and stop works just fine. """
        pass


    @pytest.mark.skip("Not implemented")
    def test_both_stop_modes_is_prohibited(self, get_pipe_manager, spec_type):
        """ Provision of both prospective and retrospective stop is bad. """
        pass


    @pytest.mark.skip("Not implemented")
    def test_complementary_specification_modes(self, get_pipe_manager):
        """ Command-line and keyword specifications can harmonize. """
        pass


    @pytest.mark.skip("Not implemented")
    def test_command_line_beats_constructor_keyword(self, get_pipe_manager):
        """ Command-line specification is favored over constructor keyword. """
        pass
