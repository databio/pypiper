""" Create dummy PipelineManager and NGSTk instance for interactive session. """

import os
from pypiper import PipelineManager
from pypiper import NGSTk


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



pm = PipelineManager(name="interactive", outfolder=os.path.expanduser("~"))
tk = NGSTk(pm=pm)
