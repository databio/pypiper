"""Pipeline base class"""

import abc
import glob
import os
from collections.abc import Callable, Iterable, Mapping
from typing import Any

from .exceptions import (
    IllegalPipelineDefinitionError,
    IllegalPipelineExecutionError,
    UnknownPipelineStageError,
)
from .manager import PipelineManager
from .stage import Stage
from .utils import (
    _checkpoint_filepath,
    _flag_name,
    _parse_stage_name,
    _translate_stage_name,
)

__all__ = ["Pipeline"]


class Pipeline(object):
    """Base class for defining a multi-stage pipeline with checkpointing.

    Example:
        class MyPipeline(Pipeline):
            def stages(self):
                return [("trim", self.trim), ("align", self.align)]

            def trim(self):
                self.manager.run("trimmomatic ...", target="trimmed.fq")

            def align(self):
                self.manager.run("bowtie2 ...", target="aligned.bam")

        p = MyPipeline(name="rnaseq", outfolder="output/")
        p.run(start_point="trim", stop_after="align")

    Requires either a PipelineManager or (name + outfolder) to construct.
    Subclasses must implement stages() returning an ordered collection.

    Args:
        name: Pipeline name. Required if manager is not provided.
        manager: Existing PipelineManager to use.
        outfolder: Output directory. Required if manager is not provided.
        args: argparse.Namespace for PipelineManager construction.
        **pl_mgr_kwargs: Additional PipelineManager keyword arguments.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(
        self,
        name: str | None = None,
        manager: PipelineManager | None = None,
        outfolder: str | None = None,
        args: Any = None,
        **pl_mgr_kwargs: Any,
    ) -> None:
        super(Pipeline, self).__init__()
        try:
            self.name = name or manager.name
        except AttributeError:
            raise TypeError(
                "If a pipeline manager isn't provided to create {}, a name is required.".format(
                    Pipeline.__name__
                )
            )
        else:
            if not self.name:
                raise ValueError(
                    "Invalid name, possible inferred from pipeline manager: {} ({})".format(
                        self.name, type(self.name)
                    )
                )

        # Determine the PipelineManager.
        if manager:
            self.manager = manager
            if outfolder:
                print(
                    "Ignoring explicit output folder ({}) and using that of "
                    "pipeline manager ({})".format(outfolder, manager.outfolder)
                )
            if name and name != manager.name:
                print(
                    "Warning: name for pipeline ('{}') doesn't match that "
                    "of the given manager ('{}')".format(name, manager.name)
                )
        elif outfolder:
            # We're guaranteed by the upfront exception block around
            # name setting that we'll either have set the name for this
            # instance to a non-null if we reach this point, and thus we're
            # protected from passing a null name argument to the pipeline
            # manager's constructor.
            self.manager = PipelineManager(self.name, outfolder, args=args, **pl_mgr_kwargs)
        else:
            raise TypeError(
                "To create a {} instance, 'manager' or 'outfolder' is required".format(
                    self.__class__.__name__
                )
            )

        # Require that checkpoints be overwritten.
        self.manager.overwrite_checkpoints = True

        # Translate stage names; do this here to complicate a hypothetical
        # attempt to override or otherwise redefine the way in which
        # stage names are handled, parsed, and translated.
        self._unordered = _is_unordered(self.stages())
        if self._unordered:
            print("NOTICE: Unordered definition of stages for pipeline {}".format(self.name))

        # Get to a sequence of pairs of key (possibly in need of translation)
        # and actual callable. Key is stage name and value is either stage
        # callable or an already-made stage object.
        stages = self.stages().items() if isinstance(self.stages(), Mapping) else self.stages()
        # Stage spec. parser handles callable validation.
        name_stage_pairs = [_parse_stage_spec(s) for s in stages]

        # Pipeline must have non-empty definition of stages.
        if not stages:
            raise IllegalPipelineDefinitionError("Empty stages")

        # Ensure that each pipeline stage is callable, and map names
        # between from external specification and internal representation.
        # We don't need to store the internal-to-external mapping, as each
        # stage will store its own name that is equivalent to the "external"
        # one, and we can use the checkpoint name derivation functions
        # to determine checkpoint name/path from stage/stage name.
        # We just use this internal-to-external mapping here, ephemerally,
        # to pretest whether there'd be a checkpoint name resolution collision.
        _internal_to_external: dict[str, str] = dict()
        self._external_to_internal: dict[str, str] = dict()
        self._stages: list[Stage] = []

        for name, stage in name_stage_pairs:
            # Use external translator to further confound redefinition.
            internal_name = _translate_stage_name(name)

            # Check that there's not a checkpoint name collision.
            if internal_name in _internal_to_external:
                already_mapped = _internal_to_external[internal_name]
                errmsg = (
                    "Duplicate stage name resolution (stage names are too "
                    "similar.) '{}' and '{}' both resolve to '{}'".format(
                        name, already_mapped, internal_name
                    )
                )
                raise IllegalPipelineDefinitionError(errmsg)

            # Store the stage name translations and the stage itself.
            self._external_to_internal[name] = internal_name
            _internal_to_external[internal_name] = name
            self._stages.append(stage)

        self.skipped, self.executed = None, None

    @property
    def outfolder(self) -> str:
        """Path to the pipeline's output folder."""
        return self.manager.outfolder

    @abc.abstractmethod
    def stages(self) -> Mapping | list:
        """Define the pipeline's processing stages. Must be overridden.

        Returns:
            Ordered collection of stages: list of (name, callable) tuples,
            list of Stage objects, or OrderedDict mapping names to callables.
        """
        pass

    @property
    def stage_names(self) -> list[str]:
        """List of stage names as defined by the pipeline author."""
        return [_parse_stage_name(s) for s in self._stages]

    def checkpoint(self, stage: Stage, msg: str = "") -> bool:
        """Write a checkpoint file for the given stage.

        Args:
            stage: Stage object to checkpoint.
            msg: Optional message for the timestamp.

        Returns:
            Whether a checkpoint file was written.
        """
        # Canonical usage model for Pipeline checkpointing through
        # implementations of this class is by automatically creating a
        # checkpoint when a conceptual unit or group of operations of a
        # pipeline completes, so fix the 'finished' parameter to the manager's
        # timestamp method to be True.
        return self.manager.timestamp(message=msg, checkpoint=stage.checkpoint_name, finished=True)

    def completed_stage(self, stage: Stage) -> bool:
        """Check whether a stage has a checkpoint file indicating completion.

        Args:
            stage: Stage to check.

        Returns:
            True if the checkpoint file exists.
        """
        check_path = _checkpoint_filepath(stage, self.manager)
        return os.path.exists(check_path)

    def halt(self, **kwargs: Any) -> None:
        """Halt the pipeline via the underlying PipelineManager."""
        self.manager.halt(**kwargs)

    def list_flags(self, only_name: bool = False) -> list[str]:
        """List flag files associated with this pipeline.

        Args:
            only_name: Return only filenames (True) or full paths (False).

        Returns:
            List of flag file paths or names.
        """
        paths = glob.glob(os.path.join(self.outfolder, _flag_name("*")))
        if only_name:
            return [os.path.split(p)[1] for p in paths]
        else:
            return paths

    def run(
        self,
        start_point: str | None = None,
        stop_before: str | None = None,
        stop_after: str | None = None,
    ) -> None:
        """Execute pipeline stages, optionally specifying start and stop points.

        Example:
            p.run()                           # run all stages
            p.run(start_point="align")        # skip to alignment
            p.run(stop_after="trim")          # stop after trimming

        Stages with existing checkpoint files from previous runs are skipped.
        When execution resumes (e.g. after a crash), stages before start_point
        are skipped entirely, and stages between start_point and the first
        uncheckpointed stage are also skipped. Cannot specify both stop_before
        and stop_after.

        Args:
            start_point: Stage name to begin execution at. Stages before this
                are skipped without checking checkpoint files.
            stop_before: Stage name to stop before (exclusive -- this stage
                does NOT run).
            stop_after: Stage name to stop after (inclusive -- this stage
                DOES run, then pipeline halts).
        """

        # Start the run with a clean slate of Stage status/label tracking.
        self._reset()

        # TODO: validate starting point against checkpoint flags for
        # TODO (cont.): earlier stages if the pipeline defines its stages as a
        # TODO (cont.): sequence (i.e., probably prohibit start point with
        # TODO (cont): nonexistent earlier checkpoint flag(s).)

        if stop_before and stop_after:
            raise IllegalPipelineExecutionError(
                "Cannot specify both inclusive and exclusive stops."
            )

        if stop_before:
            stop = stop_before
            inclusive_stop = False
        elif stop_after:
            stop = stop_after
            inclusive_stop = True
        else:
            stop = None
            inclusive_stop = None

        # Ensure that a stage name--if specified--is supported.
        for s in [start_point, stop]:
            if s is None:
                continue
            name = _parse_stage_name(s)
            if name not in self.stage_names:
                raise UnknownPipelineStageError(name, self)

        # Permit order-agnostic pipelines, but warn.
        if self._unordered and (start_point or stop_before or stop_after):
            print(
                "WARNING: Starting and stopping points are nonsense for "
                "pipeline with unordered stages."
            )

        # TODO: consider context manager based on start/stop points.

        # Determine where to start (but perhaps skip further based on
        # checkpoint completions.)
        start_index = self._start_index(start_point)
        stop_index = self._stop_index(stop, inclusive=inclusive_stop)
        assert stop_index <= len(self._stages)
        if start_index >= stop_index:
            raise IllegalPipelineExecutionError("Cannot start pipeline at or after stopping point")

        # TODO: consider storing just stage name rather than entire stage.
        # TODO (cont.): the bad case for whole-Stage is if associated data
        # TODO (cont.): (i.e., one or more args) are large.
        self.skipped.extend(self._stages[:start_index])

        # TODO: support both courses of action for non-continuous checkpoints.
        # TODO (cont.): That is, what if there's a stage with a checkpoint
        # TODO (cont.): file downstream of one without it? Naively, we'll
        # TODO (cont.): skip it, but we may want to re-run.
        skip_mode = True

        for stage in self._stages[start_index:stop_index]:
            # TODO: Note that there's no way to tell whether a non-checkpointed
            # TODO (cont.) Stage has been completed, and thus this seek
            # TODO (cont.) operation will find the first Stage, starting
            # TODO (cont.) the specified start point, either uncheckpointed or
            # TODO (cont.) for which the checkpoint file does not exist.
            # Look for checkpoint file.
            if skip_mode and self.completed_stage(stage):
                print("Skipping completed checkpoint stage: {}".format(stage))
                self.skipped.append(stage)
                continue

            # Once we've found where to being execution, ignore checkpoint
            # flags downstream if they exist since there may be dependence
            # between results from different stages.
            skip_mode = False

            print(f"Running stage: {getattr(stage, 'name', str(stage))}")

            try:
                stage.run()
            except Exception as e:
                self.manager._triage_error(e, nofail=stage.nofail)
            else:
                self.checkpoint(stage)
            self.executed.append(stage)

        # Add any unused stages to the collection of skips.
        self.skipped.extend(self._stages[stop_index:])

        # Where we stopped determines the shutdown mode.
        if stop_index == len(self._stages):
            self.wrapup()
        else:
            self.halt(raise_error=False)

    def wrapup(self) -> None:
        """Final mock stage to run after final one finishes."""
        self.manager.complete()

    def _reset(self) -> None:
        """Scrub decks with respect to Stage status/label tracking."""
        self.skipped, self.executed = [], []

    def _start_index(self, start: str | None = None) -> int:
        """Seek to the first stage to run."""
        if start is None:
            return 0
        start_stage = _translate_stage_name(start)
        internal_names = [_translate_stage_name(s.name) for s in self._stages]
        try:
            return internal_names.index(start_stage)
        except ValueError:
            raise UnknownPipelineStageError(start, self)

    def _stop_index(self, stop_point: str | None, inclusive: bool | None) -> int:
        """Determine index of stage of stopping point for run().

        Args:
            stop_point: Stopping point itself or name of it.
            inclusive: Whether the stopping point is to be regarded as
                inclusive (i.e., whether it's the final stage to run, or the one
                just beyond).

        Returns:
            Index into sequence of Pipeline's stages that indicates where to stop;
            critically, the value of the inclusive parameter here is used to
            contextualize this index such that it's always returned as an exclusive
            stopping index (i.e., execute up to the stage indexed by the value
            returned from this function).
        """
        if not stop_point:
            # Null case, no stopping point
            return len(self._stages)
        stop_name = _parse_stage_name(stop_point)
        try:
            stop_index = self.stage_names.index(stop_name)
        except ValueError:
            raise UnknownPipelineStageError(stop_name, self)
        return stop_index + 1 if inclusive else stop_index


def _is_unordered(collection: Iterable) -> bool:
    """Determine whether a collection appears to be unordered.

    This is a conservative implementation, allowing for the possibility that
    someone's implemented Mapping or Set, for example, and provided an
    __iter__ implementation that defines a consistent ordering of the
    collection's elements.

    Args:
        collection: Object to check as an unordered collection.

    Returns:
        Whether the given object appears to be unordered.

    Raises:
        TypeError: If the given "collection" is non-iterable, it's
            illogical to investigate whether it's ordered.
    """
    if not isinstance(collection, Iterable):
        raise TypeError("Non-iterable alleged collection: {}".format(type(collection)))

    return isinstance(collection, set) or isinstance(collection, dict)


def _parse_stage_spec(stage_spec: Any) -> tuple[str, Stage]:
    """Handle alternate Stage specifications, returning name and Stage.

    Isolate this parsing logic from any iteration. TypeError as single
    exception type funnel also provides a more uniform way for callers to
    handle specification errors (e.g., skip a stage, warn, re-raise, etc.)

    Args:
        stage_spec: Name and Stage tuple, or a callable.

    Returns:
        Pair of name and Stage instance from parsing input specification.

    Raises:
        TypeError: If the specification of the stage is not a supported type.
    """

    # The logic used here, a message to a user about how to specify Stage.
    req_msg = (
        "Stage specification must be either a {0} itself, a "
        "(<name>, {0}) pair, or a callable with a __name__ attribute "
        "(e.g., a non-anonymous function)".format(Stage.__name__)
    )

    # Simplest case is stage itself.
    if isinstance(stage_spec, Stage):
        return stage_spec.name, stage_spec

    # Handle alternate forms of specification.
    try:
        # Unpack pair of name and stage, requiring name first.
        name, stage = stage_spec
    except (TypeError, ValueError):
        # Normally, this sort of unpacking issue create a ValueError. Here,
        # though, we also need to catch TypeError since that's what arises
        # if an attempt is made to unpack a single function.
        # Attempt to parse stage_spec as a single named callable.
        try:
            name = stage_spec.__name__
        except AttributeError:
            raise TypeError(req_msg)
        else:
            # Control flow here indicates an anonymous function that was not
            # paired with a name. Prohibit that.
            if name == (lambda: None).__name__:
                raise TypeError(req_msg)
        stage = stage_spec

    # Ensure that the stage is callable.
    if not hasattr(stage, "__call__"):
        raise TypeError(req_msg)

    return name, Stage(stage, name=name)
