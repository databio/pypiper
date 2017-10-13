""" Pipeline base class """

import abc
from collections import Counter, OrderedDict
import sys
if sys.version_info < (3, 3):
    from collections import Sequence
else:
    from collections.abc import Sequence

from stage import translate_stage_name


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


__all__ = ["Pipeline", "UnknownPipelineStageError"]



class Pipeline(object):
    """
    Generic pipeline framework.

    :param name: Name for the pipeline; arbitrary, just used for messaging.
    :type name: str
    :param manager: The pipeline manager to
    :type manager: pypiper.PipelineManager
    :param ask_complete: How to determine whether a stage has been completed.
    :type ask_complete: function
    :param flag_complete: How to indicate to indicate to the pipeline manager
        (e.g., via the filesystem) that a stage has been completed.
    :type flag_complete: function
    """
    
    __metaclass__ = abc.ABCMeta
    
    
    def __init__(self, name, manager=None,
                 ask_complete=None, flag_complete=None):

        super(Pipeline, self).__init__()

        self.name = name
        
        if manager is None:
            if not (hasattr(ask_complete, "__call__") and 
                    hasattr(flag_complete, "__call__")):
                raise ValueError(
                    "If no pipeline manager is provided, callable stage query "
                    "and stage writing strategies must be ")
        else:
            print("Pipeline manager provided: {}".format(manager.name))
        
        self.manager = manager

        # Translate stage names; do this here to complicate a hypothetical
        # attempt to override or otherwise redefine the way in which
        # stage names are handled, parsed, and translated.
        try:
            stages = self.stages.items()
        except AttributeError:
            stages = self.stages

        if not stages:
            raise ValueError("Empty stages")

        # Check that each element is a pair, and standardize to list-of-pairs
        # representation of stage name and (hopefully callable) stage.
        try:
            name_stage_pairs = [(name, stage) for name, stage in stages]
        except ValueError:
            print(
                "Error unpacking stage name and stage itself for {}; {} stages "
                "must be defined as a name-to-callable mapping or as a "
                "collection of such pairs".format(
                    self.name, self.__class__.__name__, ))
            raise

        # Enforce stage name uniqueness.
        stage_names, _ = zip(*name_stage_pairs)    # Need non-emptiness
        stage_name_counts = Counter(stage_names)
        repeated = [(s, n) for s, n in stage_name_counts.items() if n > 1]
        if repeated:
            raise ValueError("Duplicate stage name(s): {}".format(
                ", ".join(["{} ({})".format(s, n) for s, n in repeated])))

        # Ensure that each pipeline stage is callable, and map names
        # between external specification and internal representation.
        self._internal_to_external = OrderedDict()
        self._external_to_internal = OrderedDict()
        self._stages = []
        for name, stage in name_stage_pairs:
            if not hasattr(stage, "__call__"):
                raise TypeError(
                    "Uncallable stage for stage name: {}".format(name))
            # Use external translator to further confound redefinition.
            internal_name = translate_stage_name(name)
            self._external_to_internal[name] = internal_name
            self._internal_to_external[internal_name] = name
            self._stages.append(stage)



    @abc.abstractproperty
    def stages(self):
        """
        Define the names of pipeline processing stages.

        :return: Collection of pipeline stage names.
        :rtype: Iterable[str]
        """
        pass


    @property
    def stage_names(self):
        """
        Fetch the pipeline's stage names as specified by the pipeline
        class author (i.e., not necessarily those that are used for the
        checkpoint files)

        :return:
        """
        return list(self._external_to_internal.keys())


    def run(self, start=None, stop=None):
        """
        Run the pipeline, optionally specifying start and/or stop points.

        :param start: Name of stage at which to begin execution.
        :type start: str
        :param stop: Name of stage at which to cease execution.
        :type stop: str
        """

        # TODO: validate starting point against checkpoint flags for
        # TODO (cont.): earlier stages if the pipeline defines its stages as a
        # TODO (cont.): sequence (i.e., probably prohibit start point with
        # TODO (cont): nonexistent earlier checkpoint flag(s).)

        # Ensure that a stage name--if specified--is supported.
        for s in [start, stop]:
            if not (s is None or s in self.stage_names):
                raise UnknownPipelineStageError(s, self)

        # Permit order-agnostic pipelines, but warn.
        if not isinstance(self.stages, Sequence):
            print("NOTICE: Unordered definition of stages for pipeline {}".
                  format(self.name))
            if start or stop:
                print("WARNING: Starting and stopping points are nonsense for "
                      "pipeline with unordered stages.")

        # TODO: consider context manager based on start/stop points.

        name_stage_pairs = zip(self._internal_to_external.keys(), self._stages)
        for stage_name, stage in name_stage_pairs:
            # TODO: check against start point name and for checkpoints.
            pass


    @staticmethod
    def _exec_stage(func, *args, **kwargs):
        func(*args, **kwargs)


    def is_complete(self, stage):
        """
        Determine whether the pipeline's completed the stage indicated.
        
        :param stage: Name of stage to check for completion status.
        :type stage: str
        :return: Whether this pipeline's completed the indicated stage.
        :rtype: bool
        :raises UnknownStageException: If the stage name given is undefined 
            for the pipeline, a ValueError arises.
        """
        if stage not in self.stages:
            raise UnknownPipelineStageError(stage, self)



class UnknownPipelineStageError(Exception):
    """
    Triggered by use of unknown/undefined name for a pipeline stage.
    
    :param stage_name: Name of the stage triggering the exception.
    :type stage_name: str
    :param pipeline: Pipeline for which the stage is unknown/undefined.
    :type pipeline: Pipeline
    """
    
    def __init__(self, stage_name, pipeline=None):
        message = stage_name
        if pipeline is not None:
            try:
                stages = pipeline.stages
            except AttributeError:
                # Just don't contextualize the error with known stages.
                pass
            else:
                message = "{}; defined stages: {}".\
                        format(message, ", ".join(stages))
        super(UnknownPipelineStageError, self).__init__(message)
