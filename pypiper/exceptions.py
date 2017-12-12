""" Custom pypiper exceptions """

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = ["PipelineError", "UnsupportedFiletypeException"]



class PipelineError(Exception):
    """ General pipeline error. """
    pass



class PipelineHalt(Exception):
    """
    Execution-stopping exception for halting a pipeline.

    This is useful for stopping execution of a truly script-like pipeline.
    That is, a pipeline that doesn't bundle/define stages or wrap run() calls
    in functions. In this case, we want to be able to stop the Python process
    as it chugs through a pipeline script, and we can do that by having a
    PipelineManager's halt method raise this exception.

    """
    def __init__(self, checkpoint=None, finished=None):
        if checkpoint is None:
            super(PipelineHalt, self).__init__()
        else:
            if isinstance(checkpoint, str):
                last_stage_done = checkpoint
            else:
                last_stage_done = getattr(checkpoint, "name", None) or \
                                  getattr(checkpoint, "__name__", None)
            if not last_stage_done:
                super(PipelineHalt, self).__init__()
            else:
                if finished is None:
                    msg = last_stage_done
                elif finished:
                    msg = "Finished '{}'".format(last_stage_done)
                else:
                    msg = "Stopped at '{}'".format(last_stage_done)
                super(PipelineHalt, self).__init__(msg)



class UnsupportedFiletypeException(Exception):
    """ Restrict filetype domain. """
    # Use superclass ctor to allow file name/path or extension to pass
    # through as the message for why this error is occurring.
    pass
