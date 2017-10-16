""" Conceptualize a pipeline processing phase/stage. """

import copy
from utils import pipeline_filepath

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = ["Stage"]


CHECKPOINT_EXTENSION = ".checkpoint"



class Stage(object):
    """ Single stage/phase of a pipeline; a logical processing "unit". """


    def __init__(self, func, f_args=None, f_kwargs=None,
                 name=None, checkpoint=True):
        """
        A function, perhaps with arguments, defines the stage.

        :param callable func: The processing logic that defines the stage
        :param tuple f_args: Positional arguments for func
        :param dict f_kwargs: Keyword arguments for func
        :param str name: name for the phase/stage
        """
        self.f = func
        self.f_args = f_args or tuple()
        self.f_kwargs = f_kwargs or dict()
        self.name = name or func.__name__
        self.checkpoint = checkpoint


    @property
    def checkpoint_name(self):
        """
        Determine the checkpoint name for this Stage.

        :return str | NoneType: Checkpoint name for this stage; null if this
            Stage is designated as a non-checkpoint.
        """
        return translate_stage_name(self.name) if self.checkpoint else None


    def run(self, *args, **kwargs):
        self(*args, **kwargs)


    def __call__(self, *args, **update_kwargs):
        kwargs = copy.deepcopy(self.f_kwargs)
        kwargs.update(update_kwargs)
        args = args or self.f_args
        self.f(*args, **kwargs)



def checkpoint_filename(checkpoint):
    """
    Translate a checkpoint to a filename.

    This not only adds the checkpoint file extension but also standardizes the
    way in which checkpoint names are mapped to filenames.

    :param str | Stage checkpoint: name of a pipeline phase/stage
    :return str | NoneType: standardized checkpoint name for file, plus
        extension; null if the input is a Stage that's designated as a
        non-checkpoint
    """
    if isinstance(checkpoint, Stage):
        return checkpoint.checkpoint_name
    name = translate_stage_name(checkpoint)
    return name + CHECKPOINT_EXTENSION


def checkpoint_filepath(checkpoint, pm):
    """
    Create filepath for indicated checkpoint.

    :param str checkpoint: name of a pipeline phase/stage
    :param pypiper.PipelineManager pm: manager of a pipeline instance,
        relevant here for output folder path.
    :return str: standardized checkpoint name for file, plus extension
    """
    filename = checkpoint_filename(checkpoint)
    return pipeline_filepath(pm, filename)


def translate_stage_name(stage_name):
    """
    Account for potential variability in stage/phase name definition.

    Since a pipeline author is free to name his/her processing phases/stages
    as desired, but these choices influence file names, enforce some
    standardization. Specifically, prohibit potentially problematic spaces.

    :param stage_name: Name of the pipeline phase/stage.
    :type stage_name: str
    :return: Standardized pipeline phase/stage name.
    :rtype: str
    """
    # Cast to string to ensure that indexed stages (ints are handled).
    return str(stage_name).lower().replace(" ", "-")

