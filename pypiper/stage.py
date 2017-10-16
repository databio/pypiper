""" Conceptualize a pipeline processing phase/stage. """

import copy
from utils import translate_stage_name

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

        :return: Checkpoint name for this stage; null if this Stage is
            designated as a non-checkpoint.
        :rtype: str | NoneType
        """
        return translate_stage_name(self.name) if self.checkpoint else None


    def run(self, *args, **kwargs):
        """ Alternate form for direct call; execute stage. """
        self(*args, **kwargs)


    def __call__(self, *args, **update_kwargs):
        """ Execute the stage, allowing updates to args/kwargs. """
        kwargs = copy.deepcopy(self.f_kwargs)
        kwargs.update(update_kwargs)
        args = args or self.f_args
        self.f(*args, **kwargs)


    def __repr__(self):
        return "{klass} '{n}': f={f}, args={pos}, kwargs={kwd}, " \
               "checkpoint={check}".format(klass=self.__class__.__name__,
                f=self.f, n=self.name, pos=self.f_args, kwd=self.f_kwargs,
                check=self.checkpoint)


    def __str__(self):
        return "{}: '{}'".format(self.__class__.__name__, self.name)



def checkpoint_filename(checkpoint):
    """
    Translate a checkpoint to a filename.

    This not only adds the checkpoint file extension but also standardizes the
    way in which checkpoint names are mapped to filenames.

    :param checkpoint: name of a pipeline phase/stage
    :type checkpoint: str | Stage
    :return: standardized checkpoint name for file, plus extension;
        null if the input is a Stage that's designated as a non-checkpoint
    :rtype: str | NoneType
    """
    if isinstance(checkpoint, Stage):
        base = checkpoint.checkpoint_name
    else:
        base = translate_stage_name(checkpoint)
    return base + CHECKPOINT_EXTENSION
