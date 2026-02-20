"""Custom pypiper exceptions"""

from typing import Any

__all__ = [
    "PipelineError",
    "PipelineHalt",
    "IllegalPipelineDefinitionError",
    "IllegalPipelineExecutionError",
    "MissingCheckpointError",
    "UnknownPipelineStageError",
    "UnsupportedFiletypeException",
    "SubprocessError",
]


class PipelineError(Exception):
    """General pipeline error."""

    pass


class SubprocessError(Exception):
    """Error from a failed subprocess call."""

    def __init__(self, msg="", returncode=None, cmd=None):
        self.returncode = returncode
        self.cmd = cmd
        super().__init__(msg)


class IllegalPipelineDefinitionError(PipelineError):
    pass


class IllegalPipelineExecutionError(PipelineError):
    """Represent cases of illogical start/stop run() declarations."""

    pass


class MissingCheckpointError(Exception):
    """Represent case of expected but absent checkpoint file."""

    def __init__(self, checkpoint: str, filepath: str) -> None:
        msg = (
            "Checkpoint '{checkpoint}' not found at expected path: '{filepath}'. "
            "This checkpoint file is created when the '{checkpoint}' stage completes successfully. "
            "The stage may not have been run yet, or its output was deleted. "
            "To re-run from this stage, use start_point='{checkpoint}' or re-run the full pipeline "
            "with new_start=True (CLI: -N).".format(checkpoint=checkpoint, filepath=filepath)
        )
        super(MissingCheckpointError, self).__init__(msg)


class UnknownPipelineStageError(Exception):
    """Triggered by use of unknown/undefined name for a pipeline stage.

    Args:
        stage_name: Name of the stage triggering the exception.
        pipeline: Pipeline for which the stage is unknown/undefined.
    """

    def __init__(self, stage_name: str, pipeline: Any = None) -> None:
        message = "Unknown pipeline stage: '{stage}'.".format(stage=stage_name)
        if pipeline is not None:
            try:
                stages = pipeline.stages()
            except AttributeError:
                # Just don't contextualize the error with known stages.
                pass
            else:
                stage_names = ", ".join(map(str, stages))
                message = (
                    "Unknown pipeline stage: '{stage}'. "
                    "Available stages are: [{stage_names}]. "
                    "Check for typos in the stage name passed to start_point, stop_before, or stop_after.".format(
                        stage=stage_name, stage_names=stage_names
                    )
                )
        super(UnknownPipelineStageError, self).__init__(message)


class PipelineHalt(Exception):
    """
    Execution-stopping exception for halting a pipeline.

    This is useful for stopping execution of a truly script-like pipeline.
    That is, a pipeline that doesn't bundle/define stages or wrap run() calls
    in functions. In this case, we want to be able to stop the Python process
    as it chugs through a pipeline script, and we can do that by having a
    PipelineManager's halt method raise this exception.

    """

    def __init__(self, checkpoint: str | None = None, finished: bool | None = None) -> None:
        if checkpoint is None:
            super(PipelineHalt, self).__init__()
        else:
            if isinstance(checkpoint, str):
                last_stage_done = checkpoint
            else:
                last_stage_done = getattr(checkpoint, "name", None) or getattr(
                    checkpoint, "__name__", None
                )
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
    """Restrict filetype domain."""

    # Use superclass ctor to allow file name/path or extension to pass
    # through as the message for why this error is occurring.
    pass
