""" Conceptualize a pipeline processing phase/stage. """

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

        :param func:
        :param f_args:
        :param f_kwargs:
        :param str name: name for the phase/stage
        """
        self.f = func
        self.f_args = f_args or tuple()
        self.f_kwargs = f_kwargs or dict()
        # Currently unused, but a hook for defining checkpoint filename
        # different than
        self.name = name
        self.checkpoint = checkpoint


    def __call__(self, *args, **kwargs):
        self.f(*self.f_args, **self.f_kwargs)



def checkpoint_filename(checkpoint):
    """
    Translate a checkpoint to a filename.

    This not only adds the checkpoint file extension but also standardizes the
    way in which checkpoint names are mapped to filenames.

    :param str checkpoint: name of a pipeline phase/stage
    :return str: standardized checkpoint name for file, plus extension
    """
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
