""" Tests for the Pipeline data type """

from functools import partial
import glob
import os
import pytest
from pypiper import Pipeline
from pypiper.manager import COMPLETE_FLAG, PAUSE_FLAG, RUN_FLAG
from pypiper.pipeline import \
        checkpoint_filepath, IllegalPipelineDefinitionError, \
        IllegalPipelineExecutionError, UnknownPipelineStageError
from pypiper.stage import Stage
from pypiper.utils import \
    flag_name, pipeline_filepath, checkpoint_filename, translate_stage_name
from .helpers import named_param
from tests.conftest import \
    write_file1, write_file2, write_file3, \
    CONTENTS, FILENAMES, FILE1_NAME, FILE_TEXT_PAIRS, \
    OUTPUT_SUFFIX, TEST_PIPE_NAME


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



BASIC_ACTIONS = [write_file1, write_file2, write_file3]
STAGE_SPECS = ["stage", "name", "function"]



def pytest_generate_tests(metafunc):
    """ Dynamic creation and parameterization of cases in this module. """
    if "test_type" in metafunc.fixturenames and \
            metafunc.cls == MostBasicPipelineTests:
        metafunc.parametrize(
                argnames="test_type", 
                argvalues=["effects", "stage_labels", 
                           "checkpoints", "pipe_flag"])



@pytest.fixture
def stage(request):
    """
    Provide a test case with a pipeline stage.

    :param pytest.fixtures.FixtureRequest request: Test case in need of the
        parameterization
    :return str | pypiper.Stage | function: The stage provided, in the
        desired format.
    """
    # Ensure that the test case has the needed fixtures.
    stage_fixture_name = "point"
    spec_type_name = "spec_type"
    for fixture_name in [stage_fixture_name, spec_type_name]:
        if fixture_name not in request.fixturenames:
            raise ValueError("No '{}' fixture".format(stage_fixture_name))
    s = request.getfixturevalue(stage_fixture_name)
    spec_type = request.getfixturevalue(spec_type_name)
    return _parse_stage(s, spec_type)



def _parse_stage(s, spec_type):
    """
    Do a type transformation on a Stage function.

    :param callable s: the callable hypothetically backing a Stage
    :param str spec_type: indication of the type transformation to perform
    :return object: representation of Stage, of requested type
    """
    if spec_type == "stage":
        s = Stage(s)
    elif spec_type == "name":
        s = Stage(s).name
    elif spec_type == "function":
        assert hasattr(s, "__call__")
    else:
        raise ValueError("Unknown specification type: {}".format(spec_type))
    return s



class EmptyStagesPipeline(Pipeline):
    """ Illegal (via empty stages) Pipeline definition. """

    def __init__(self, manager):
        super(EmptyStagesPipeline, self).__init__(
                TEST_PIPE_NAME, manager=manager)

    def stages(self):
        return []


class NameCollisionPipeline(Pipeline):
    """ Illegal (via empty stages) Pipeline definition. """

    def __init__(self, manager):
        super(NameCollisionPipeline, self).__init__(
                TEST_PIPE_NAME, manager=manager)

    def stages(self):
        name = "write file1"
        return [("write file1", write_file1),
                (translate_stage_name(name), write_file1)]



class RunPipelineCornerCaseTests:
    """ Tests for exceptional cases of pipeline execution. """


    @named_param(argnames="point", argvalues=BASIC_ACTIONS)
    @named_param(argnames="spec_type", argvalues=STAGE_SPECS)
    @named_param(argnames="inclusive", argvalues=[False, True])
    def test_start_point_equals_stop(
            self, dummy_pipe, point, spec_type, stage, inclusive):
        """ Start=stop is only permitted if stop should be run. """

        _assert_pipeline_initialization(dummy_pipe)

        # Inclusion determines how to make the call, and the expectation.
        if inclusive:
            # start_point = inclusive stop --> single stage runs.
            dummy_pipe.run(start_point=stage, stop_after=stage)
            _assert_checkpoints(dummy_pipe, [stage])
        else:
            # start_point = exclusive stop --> exception
            with pytest.raises(IllegalPipelineExecutionError):
                dummy_pipe.run(start_point=stage, stop_before=stage)


    @pytest.mark.parametrize(
            argnames=["start_point", "stop"],
            argvalues=[(write_file2, write_file1),
                       (write_file3, write_file2),
                       (write_file3, write_file1)])
    @pytest.mark.parametrize(argnames="spec_type", argvalues=STAGE_SPECS)
    @pytest.mark.parametrize(
            argnames="stop_type", argvalues=["stop_before", "stop_after"])
    def test_start_point_after_stop(
            self, dummy_pipe, start_point, stop, stop_type, spec_type):
        """ Regardless of specification type, start > stop is prohibited. """
        start_point = _parse_stage(start_point, spec_type)
        stop = _parse_stage(stop, spec_type)
        with pytest.raises(IllegalPipelineExecutionError):
            dummy_pipe.run(**{"start_point": start_point, stop_type: stop})


    @named_param(
            argnames="undefined_stage",
            argvalues=["unsupported-pipeline-stage", "unknown_phase"])
    @named_param(argnames="stage_point",
                 argvalues=["start_point", "stop_before", "stop_after"])
    def test_unknown_stage(self, dummy_pipe, undefined_stage, stage_point):
        """ Start specification must be of known stage name. """
        with pytest.raises(UnknownPipelineStageError):
            dummy_pipe.run(**{stage_point: undefined_stage})


    @named_param(argnames="stop_before", argvalues=BASIC_ACTIONS)
    @named_param(argnames="stop_after", argvalues=BASIC_ACTIONS)
    @named_param(argnames="spec_type", argvalues=STAGE_SPECS)
    def test_stop_before_and_stop_after(
            self, dummy_pipe, stop_before, stop_after, spec_type):
        """ Inclusive and exclusive stop cannot both be provided. """
        inclusive_stop = _parse_stage(stop_after, spec_type)
        exclusive_stop = _parse_stage(stop_before, spec_type)
        kwargs = {"stop_before": exclusive_stop, "stop_after": inclusive_stop}
        with pytest.raises(IllegalPipelineExecutionError):
            dummy_pipe.run(**kwargs)


    def test_empty_stages_is_prohibited(self, pl_mgr):
        """ Pipeline must have non-empty stages """
        with pytest.raises(IllegalPipelineDefinitionError):
            EmptyStagesPipeline(manager=pl_mgr)


    def test_stage_name_collision_is_prohibited(self, pl_mgr):
        """ Each stage needs unique translation, used for checkpoint file. """
        with pytest.raises(IllegalPipelineDefinitionError):
            NameCollisionPipeline(manager=pl_mgr)



class MostBasicPipelineTests:
    """ Test pipeline defined with notion of 'absolute minimum' config. """


    def test_runs_through_full(self, dummy_pipe, test_type):
        """ The entire basic pipeline should execute. """

        # Start with clean output folder.
        _assert_pipeline_initialization(dummy_pipe)

        # Make the call under test.
        dummy_pipe.run(start_point=None, stop_before=None, stop_after=None)
        
        if test_type == "effects":
            # We're interested in existence and content of targets.
            exp_files, _ = zip(*FILE_TEXT_PAIRS)
            _assert_output(dummy_pipe, exp_files)
            fpath_text_pairs = [(pipeline_filepath(dummy_pipe, fname), content)
                                for fname, content in FILE_TEXT_PAIRS]
            for fpath, content in fpath_text_pairs:
                _assert_expected_content(fpath, content)

        elif test_type == "checkpoints":
            # Interest is on checkpoint file existence.
            for stage in dummy_pipe.stages():
                chkpt_fpath = checkpoint_filepath(stage, dummy_pipe)
                try:
                    assert os.path.isfile(chkpt_fpath)
                except AssertionError:
                    print("Stage '{}' file doesn't exist: '{}'".format(
                        stage.name, chkpt_fpath))
                    raise

        elif test_type == "stage_labels":
            # Did the Pipeline correctly track it's execution?
            _assert_stage_states(dummy_pipe, [], dummy_pipe.stages())

        elif test_type == "pipe_flag":
            # The final flag should be correctly set.
            _assert_pipeline_completed(dummy_pipe)

        else:
            raise ValueError("Unknown test type: {}".format(test_type))


    def test_skip_completed(self, dummy_pipe, test_type):
        """ Pre-completed stage(s) are skipped. """

        _assert_pipeline_initialization(dummy_pipe)

        first_stage = dummy_pipe.stages()[0]
        first_stage_chkpt_fpath = checkpoint_filepath(first_stage, dummy_pipe)
        open(first_stage_chkpt_fpath, 'w').close()
        assert os.path.isfile(first_stage_chkpt_fpath)

        exp_skips = [first_stage]
        exp_execs = dummy_pipe.stages()[1:]

        # This should neither exist nor be created.
        first_stage_outfile = pipeline_filepath(
                dummy_pipe.manager, filename=FILE1_NAME)
        assert not os.path.isfile(first_stage_outfile)
        
        # Do the action.
        dummy_pipe.run()
        
        if test_type == "effects":
            # We should not have generated the first stage's output file.
            # That notion is covered in the outfiles assertion.
            _assert_output(dummy_pipe, FILENAMES[1:])
            assert not os.path.isfile(first_stage_outfile)
            # We should have the correct content in files from stages run.
            for fname, content in FILE_TEXT_PAIRS[1:]:
                fpath = os.path.join(dummy_pipe.outfolder, fname)
                _assert_expected_content(fpath, content)

        elif test_type == "checkpoints":
            # We've manually created the first checkpoint file, and the
            # others should've been created by the run() call.
            _assert_checkpoints(dummy_pipe, exp_stages=dummy_pipe.stages())

        elif test_type == "stage_labels":
            _assert_stage_states(dummy_pipe, exp_skips, exp_execs)

        elif test_type == "pipe_flag":
            _assert_pipeline_completed(dummy_pipe)
        else:
            raise ValueError("Unknown test type: '{}'".format(test_type))


    @named_param(argnames="start_index",
                 argvalues=range(len(BASIC_ACTIONS) - 1))
    @named_param(argnames="start_spec_type",
                 argvalues=["stage", "function", "name"])
    def test_execution_allows_specific_starting_point(
            self, dummy_pipe, test_type, start_index, start_spec_type):
        """ A pipeline may be started from an arbitrary checkpoint. """
        _assert_pipeline_initialization(dummy_pipe)
        s = _parse_stage(BASIC_ACTIONS[start_index], start_spec_type)
        dummy_pipe.run(start_point=s)
        if test_type == "effects":
            exp_files = FILENAMES[start_index:]
            _assert_output(dummy_pipe, exp_files)
            fpaths = [pipeline_filepath(dummy_pipe.manager, filename=fn)
                      for fn in exp_files]
            for fp, content in zip(fpaths, CONTENTS[start_index:]):
                _assert_expected_content(fp, content)
        elif test_type == "checkpoints":
            # Ensure exact collection of checkpoint files (no more, no less).
            _assert_checkpoints(dummy_pipe, BASIC_ACTIONS[start_index:])
        elif test_type == "stage_labels":
            # Ensure match between skipped and executed stage expectations
            # and observations.
            _assert_stage_states(dummy_pipe, BASIC_ACTIONS[:start_index],
                                 BASIC_ACTIONS[start_index:])
        elif test_type == "pipe_flag":
            _assert_pipeline_completed(dummy_pipe)
        else:
            raise ValueError("Unknown test type: '{}'".format(test_type))


    def test_all_checkpoints_after_first_executed_are_overwritten(
            self, dummy_pipe, test_type):
        """ Potential for dependent results means execution is contiguous. """

        # Start fresh.
        _assert_pipeline_initialization(dummy_pipe)

        # Create checkpoint files for all but first stage.
        fpath_time_pairs = []
        for s in BASIC_ACTIONS[1:]:
            check_fpath = checkpoint_filepath(s, dummy_pipe.manager)
            open(check_fpath, 'w').close()
            fpath_time_pairs.append((check_fpath, os.path.getmtime(check_fpath)))
            assert os.path.isfile(check_fpath)

        # Ensure that we have checkpoint file for each but first stage.
        _assert_checkpoints(dummy_pipe, BASIC_ACTIONS[1:])
        # We should have output for no stage.
        _assert_output(dummy_pipe, [])

        # Make the call under test.
        dummy_pipe.run()

        if test_type == "effects":
            _assert_output(dummy_pipe, FILENAMES)
        elif test_type == "checkpoints":
            _assert_checkpoints(dummy_pipe, BASIC_ACTIONS)
        elif test_type == "stage_labels":
            _assert_stage_states(dummy_pipe, expected_skipped=[],
                                 expected_executed=BASIC_ACTIONS)
        elif test_type == "pipe_flag":
            _assert_pipeline_completed(dummy_pipe)
        else:
            raise ValueError("Unknown test type: {}".format(test_type))


    @named_param(argnames="stop_index", argvalues=range(1, len(BASIC_ACTIONS)))
    @named_param(argnames="spec_type", argvalues=STAGE_SPECS)
    @named_param(argnames="stop_type", argvalues=["stop_before", "stop_after"])
    def test_stop(self, dummy_pipe, test_type, stop_index, spec_type, stop_type):
        """ A pipeline is capable of halting at/after a specified stage. """

        # Negative control / pretest.
        _assert_pipeline_initialization(dummy_pipe)

        # Get the stop point in the correct format.
        stop = _parse_stage(BASIC_ACTIONS[stop_index], spec_type)

        # Make the call under test.
        dummy_pipe.run(**{stop_type: stop})

        # For forming expectations, indexing is exclusive.
        # So if the initial specification was inclusive, we need to
        # increment our expectation-indexing bound.
        if stop_type == "stop_after":
            stop_index += 1

        if test_type == "effects":
            exp_files = FILENAMES[:stop_index]
            _assert_output(dummy_pipe, exp_files)
            fpaths = [pipeline_filepath(dummy_pipe.manager, filename=fn)
                      for fn in exp_files]
            for fp, content in zip(fpaths, CONTENTS[:stop_index]):
                _assert_expected_content(fp, content)

        elif test_type == "checkpoints":
            _assert_checkpoints(dummy_pipe, BASIC_ACTIONS[:stop_index])
        elif test_type == "stage_labels":
            _assert_stage_states(
                    dummy_pipe, expected_skipped=BASIC_ACTIONS[stop_index:],
                    expected_executed=BASIC_ACTIONS[:stop_index])
        elif test_type == "pipe_flag":
            if (stop_index == len(BASIC_ACTIONS)) and \
                    (stop_type == "stop_after"):
                _assert_pipeline_completed(dummy_pipe)
            else:
                _assert_pipeline_halted(dummy_pipe)
        else:
            raise ValueError("Unknown test type: '{}'".format(test_type))



@named_param(
    argnames="spec_type",
    argvalues=["filename", "filepath", "stage", "stage_name"])
@named_param(argnames="completed", argvalues=[False, True])
def test_stage_completion_determination(dummy_pipe, spec_type, completed):
    """ Pipeline responds to variety of request forms of checkpoint status. """

    # Allow dummy stage definition and determination of filename.
    def dummy_test_func():
        pass

    chkpt_name = checkpoint_filename(
            dummy_test_func.__name__, pipeline_name=dummy_pipe.name)
    chkpt_fpath = checkpoint_filepath(chkpt_name, dummy_pipe.manager)

    # Determine how to request the checkpoint completion status.
    if spec_type == "filename":
        s = chkpt_name
    elif spec_type == "filepath":
        s = chkpt_fpath
    elif spec_type in ["stage", "stage_name"]:
        s = Stage(dummy_test_func)
        if spec_type == "stage_name":
            s = s.name
    else:
        raise ValueError("Unknown test spec type: {}".format(spec_type))

    # Touch the checkpoint file iff we're positively testing completion.
    if completed:
        open(chkpt_fpath, 'w').close()

    # Check the completion status for concordance with expectation.
    # Print a bit of info to inform hypotheses about the source of a
    # hypothetical test error/failure.
    outfolder_contents = os.listdir(dummy_pipe.outfolder)
    print("Pipeline outfolder contents: {}".format(outfolder_contents))
    print("Contents as pipeline files: {}".format(
        [pipeline_filepath(dummy_pipe.manager, f) for f in outfolder_contents]))
    print("Checking completion status: {} ({})".format(s, type(s)))
    observed_completion = dummy_pipe.completed_stage(s)
    if completed:
        assert observed_completion
    else:
        assert not observed_completion



def _assert_checkpoints(pl, exp_stages):
    """
    Assert equivalence between expected and observed checkpoint files.

    :param pypiper.Pipeline pl: Pipeline for which to validate checkpoints.
    :param Iterable[str | pypiper.Stage] exp_stages: Stage(s) or name(s)
        that are expected to have checkpoint file present.
    """
    exp_fpaths = [checkpoint_filepath(s, pl.manager) for s in exp_stages]
    obs_fpaths = glob.glob(checkpoint_filepath("*", pm=pl.manager))
    assert len(exp_fpaths) == len(obs_fpaths)
    assert set(exp_fpaths) == set(obs_fpaths)



def _assert_expected_content(fpath, content):
    """
    Determine whether a filepath has the expected content.

    :param str fpath: Path to file of interest.
    :param str content: Expected file content.
    :return bool: Whether observation matches expectation.
    """
    assert os.path.isfile(fpath)
    exp_content = content.split(os.linesep)
    with open(fpath, 'r') as f:
        obs_content = [l.rstrip(os.linesep) for l in f.readlines()]
    assert exp_content == obs_content



def _assert_output(pl, expected_filenames):
    """
    Assert equivalence--with respect to presence only--between expected
    collection of output file and the observed output file collection for a
    pipeline.

    :param pypiper.Pipeline pl: pipeline for which output is to be checked
    :param Iterable[str] expected_filenames:
    :return:
    """
    obs_outfiles = glob.glob(pipeline_filepath(
            pl.manager, "*{}".format(OUTPUT_SUFFIX)))
    assert len(expected_filenames) == len(obs_outfiles)
    expected_filepaths = []
    for fname in expected_filenames:
        fpath = fname if os.path.isabs(fname) else \
                pipeline_filepath(pl.manager, filename=fname)
        expected_filepaths.append(fpath)
    assert set(expected_filepaths) == set(obs_outfiles)



def _assert_pipeline_status(pl, flag):
    """ Assert, based on flag file presence, that a pipeline's completed. """
    flags = glob.glob(pipeline_filepath(pl.manager, filename=flag_name("*")))
    assert 1 == len(flags)
    exp_flag = pipeline_filepath(pl, suffix="_" + flag_name(flag))
    try:
        assert os.path.isfile(exp_flag)
    except AssertionError:
        print("FLAGS: {}".format(flags))
        raise



_assert_pipeline_completed = partial(
        _assert_pipeline_status, flag=COMPLETE_FLAG)
_assert_pipeline_halted = partial(_assert_pipeline_status, flag=PAUSE_FLAG)



def _assert_pipeline_initialization(pl):
    """
    Assert that a test case begins with output folder in expected state.

    :param pypiper.Pipeline pl: Pipeline instance for test case.
    """
    # TODO: implement.
    suffices = {"_commands.sh", "_profile.tsv",
                "_{}".format(flag_name(RUN_FLAG))}
    exp_init_contents = \
            [pipeline_filepath(pl.manager, suffix=s) for s in suffices]
    obs_init_contents = [pipeline_filepath(pl.manager, filename=n)
                         for n in os.listdir(pl.outfolder)]
    assert len(exp_init_contents) == len(obs_init_contents)
    assert set(exp_init_contents) == set(obs_init_contents)



def _assert_stage_states(pl, expected_skipped, expected_executed):
    """ Assert equivalence between expected and observed stage states. """
    def _ensure_stage(s):
        return s if isinstance(s, Stage) else Stage(s)
    expected_skipped = [_ensure_stage(s) for s in expected_skipped]
    expected_executed = [_ensure_stage(s) for s in expected_executed]
    assert expected_skipped == pl.skipped
    assert expected_executed == pl.executed
