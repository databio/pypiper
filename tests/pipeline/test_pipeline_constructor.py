""" Tests for construction of a Pipeline """

import pytest

from pypiper import Pipeline, PipelineManager
from tests.helpers import named_param


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def test_pipeline_requires_no_manager():
    """ Pipeline construction doesn't require pipeline manager argument. """

    # Ability to create pipeline without an exception implies lack of
    # requirement for a manager.
    _EmptyPipeline("empty")



def test_pipeline_requires_either_manager_or_outfolder():
    """ Pipeline must be passed pipeline manager or output folder. """
    with pytest.raises(TypeError):
        _EmptyPipeline()



def test_pipeline_creates_manager():
    """ If not passed a pipeline manager, a pipeline creates one. """
    name = "empty"
    empty = _EmptyPipeline(name)
    pm = empty.manager
    assert isinstance(pm, PipelineManager)



@named_param()
def test_pipeline_uses_manager_outfolder():
    """ The pipeline's output folder is that of its manager. """
    empty =



def test_created_manager_name_matches_pipeline_name():
    pass



def test_pipeline_ignores_outfolder_if_manager_is_passed():
    """ Manager's outfolder trumps explicit outfolder if both are passed. """
    pass



def test_pipeline_requires_stages_definition():
    """ To create a pipeline, define stages (execution steps). """

    class NoStagesPipeline(Pipeline):
        pass

    name = "test-pipe"

    # Sensitivity: test exception for bad case.
    with pytest.raises(TypeError):
        NoStagesPipeline(name)
    # Specificity: test no exception for good case.
    with pytest.raises(None):
        _EmptyPipeline(name)



@pytest.fixture
def empty_pipeline(request):
    """ Provide test case with minimal pipeline instance. """
    if "pipe_name" in request.fixturenames:
        name = request.getfixturevalue("pipe_name")
    else:
        name = "empty"
    return _EmptyPipeline(name)



class _EmptyPipeline(Pipeline):
    """ Minimal pipeline declaration. """

    def stages(self):
        pass
