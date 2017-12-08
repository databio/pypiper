""" Test effects of construction of a pipeline manager. """

import argparse
import pytest
from pypiper.manager import CHECKPOINT_SPECIFICATIONS
from tests.helpers import named_param


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def pytest_generate_tests(metafunc):
    """ Dynamic test case generation for this module's test cases. """
    if "spec_type" in metafunc.fixturenames:
        metafunc.parametrize(
                argnames="spec_type", argvalues=["cmdl", "ctor"])



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


    @named_param("start_point", ["merge_input", "filter_reads"])
    @named_param("stop_point", ["align_reads", "calc_stats"])
    @named_param("stop_type", ["stop_before", "stop_after"])
    def test_start_and_stop(self, get_pipe_manager, spec_type,
                            stop_type, start_point, stop_point):
        """ Specifying both start and stop works just fine. """
        spec_data = {"start_point": start_point, stop_type: stop_point}
        if spec_type == "cmdl":
            kwargs = {"args": argparse.Namespace(**spec_data)}
        else:
            kwargs = spec_data
        pm = get_pipe_manager(name="start-and-stop-test", **kwargs)
        assert start_point == pm.start_point
        assert stop_point == getattr(pm, stop_type)


    @named_param("stop_before", ["align_reads", "call_peaks"])
    @named_param("stop_after", ["fastqc", "align_reads"])
    @named_param("stop_before_type", ["cmdl", "ctor"])
    @named_param("stop_after_type", ["cmdl", "ctor"])
    def test_both_stop_modes_is_prohibited(
            self, get_pipe_manager, stop_before_type,
            stop_after_type, stop_before, stop_after):
        """ Provision of both prospective and retrospective stop is bad. """
        raw_kwargs = {"stop_before": stop_before, "stop_after": stop_after}
        cmdl_kwargs = {}
        if stop_before_type == "cmdl":
            cmdl_kwargs["stop_before"] = raw_kwargs.pop("stop_before")
        if stop_after_type == "cmdl":
            cmdl_kwargs["stop_after"] = raw_kwargs.pop("stop_after")
        args = argparse.Namespace(**cmdl_kwargs)
        with pytest.raises(TypeError):
            get_pipe_manager(name="test-double-stop", args=args, **raw_kwargs)


    @pytest.mark.parametrize(
        argnames=["start_point", "stop_point"],
        argvalues=[("fastqc", "align_reads"), ("align_reads", "call_peaks")])
    @pytest.mark.parametrize(
        argnames=["start_spec_type", "stop_spec_type"],
        argvalues=[("cmdl", "ctor"), ("ctor", "cmdl")])
    @named_param("stop_type", ["stop_before", "stop_after"])
    def test_complementary_specification_modes(
            self, get_pipe_manager, start_spec_type, stop_spec_type,
            stop_type, start_point, stop_point):
        """ Command-line and keyword specifications can harmonize. """
        raw_kwargs = {"start_point": start_point, stop_type: stop_point}
        cmdl_kwargs = {}
        if start_spec_type == "cmdl":
            cmdl_kwargs["start_point"] = raw_kwargs.pop("start_point")
        if stop_spec_type == "cmdl":
            cmdl_kwargs[stop_type] = raw_kwargs.pop(stop_type)
        args = argparse.Namespace(**cmdl_kwargs)
        pm = get_pipe_manager(name="complementary-test",
                              args=args, **raw_kwargs)
        assert start_point == pm.start_point
        assert stop_point == getattr(pm, stop_type)


    @named_param(
        "check_specs",
        [["start_point"], ["stop_before"], ["stop_after"],
         ["start_point", "stop_before"], ["start_point", "stop_after"]])
    def test_command_line_beats_constructor_keyword(
            self, get_pipe_manager, check_specs):
        """ Command-line specification is favored over constructor keyword. """

        # Declare values to use for respective specification modes.
        cmdl_values = {"start_point": "merge_input",
                       "stop_before": "call_peaks",
                       "stop_after": "align_reads"}
        ctor_values = {"start_point": "fastqc",
                       "stop_before": "align_reads",
                       "stop_after": "filter_reads"}

        # Create specifications based on current test case parameterization.
        cmdl_kwargs ={cp_spec: cmdl_values[cp_spec]
                      for cp_spec in check_specs}
        ctor_kwargs = {cp_spec: ctor_values[cp_spec]
                       for cp_spec in check_specs}
        args = argparse.Namespace(**cmdl_kwargs)

        # Build the pipeline manager.
        pm = get_pipe_manager(name="cmdl-preference", args=args, **ctor_kwargs)

        # Verify the preference for command-line value over variable keyword
        # argument value.
        for cp_spec in check_specs:
            assert cmdl_kwargs[cp_spec] == getattr(pm, cp_spec)

        # Verify that the non-specified values were set to null.
        for cp_spec in set(CHECKPOINT_SPECIFICATIONS) - set(check_specs):
            assert getattr(pm, cp_spec) is None
