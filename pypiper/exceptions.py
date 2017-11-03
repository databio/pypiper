""" Custom pypiper exceptions """

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = ["PipelineError", "UnsupportedFiletypeException"]



class PipelineError(Exception):
    """ General pipeline error. """
    pass



class UnsupportedFiletypeException(Exception):
    """ Restrict filetype domain. """
    # Use superclass ctor to allow file name/path or extension to pass
    # through as the message for why this error is occurring.
    pass
