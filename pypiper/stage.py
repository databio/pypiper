""" Conceptualize a pipeline processing phase/stage. """

import abc
from utils import pipeline_filepath

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



CHECKPOINT_EXTENSION = ".checkpoint"



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
    return stage_name.lower().replace(" ", "-")



class Resumable(object):
    """ An executor that supports the notion of 'resuming.' """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def resume(self):
        """ Resume a pipeline. """
        pass



class Restartable(object):
    """ An executor that supports the notion of 'restarting.' """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def restart(self):
        """ Restart a pipeline """
        pass



class PipelineStage(object):
    """ A generic pipeline processing phase/stage """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        pass
