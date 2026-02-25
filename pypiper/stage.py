"""Conceptualize a pipeline processing phase/stage."""

import copy
from collections.abc import Callable
from typing import Any

from .utils import _translate_stage_name

__all__ = ["Stage"]


class Stage(object):
    """A single checkpointable unit of pipeline processing.

    Example:
        stage = Stage(my_function, name="alignment", nofail=True)
        stage()  # executes my_function

    Args:
        func: Callable that implements the stage's logic.
        f_args: Positional arguments for func.
        f_kwargs: Keyword arguments for func.
        name: Stage name. Defaults to func.__name__.
        checkpoint: Whether to create a checkpoint file. Default: True.
        nofail: Allow stage failure without failing the pipeline.
    """

    def __init__(
        self,
        func: Callable,
        f_args: tuple | None = None,
        f_kwargs: dict | None = None,
        name: str | None = None,
        checkpoint: bool = True,
        *,
        nofail: bool = False,
    ) -> None:
        """Create a Stage from a callable."""
        if isinstance(func, Stage):
            raise TypeError("Cannot create Stage from Stage")
        super(Stage, self).__init__()
        self.f = func
        self.f_args = f_args or tuple()
        self.f_kwargs = f_kwargs or dict()
        self.name = name or func.__name__
        self.checkpoint = checkpoint
        self.nofail = nofail

    @property
    def checkpoint_name(self) -> str | None:
        """Checkpoint file name for this stage, or None if not a checkpoint."""
        return _translate_stage_name(self.name) if self.checkpoint else None

    def run(self, *args: Any, **kwargs: Any) -> None:
        """Execute the stage (alias for direct call)."""
        self(*args, **kwargs)

    def __call__(self, *args: Any, **update_kwargs: Any) -> None:
        """Execute the stage, allowing updates to args/kwargs."""
        kwargs = copy.deepcopy(self.f_kwargs)
        kwargs.update(update_kwargs)
        args = args or self.f_args
        self.f(*args, **kwargs)

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Stage)
            and self.f.__name__ == other.f.__name__
            and (
                {k: v for k, v in self.__dict__.items() if k != "f"}
                == {k: v for k, v in other.__dict__.items() if k != "f"}
            )
        )

    def __ne__(self, other: object) -> bool:
        return not (self == other)

    def __repr__(self) -> str:
        return "{klass} '{n}': f={f}, args={pos}, kwargs={kwd}, checkpoint={check}".format(
            klass=self.__class__.__name__,
            f=self.f,
            n=self.name,
            pos=self.f_args,
            kwd=self.f_kwargs,
            check=self.checkpoint,
        )

    def __str__(self) -> str:
        return "{}: '{}'".format(self.__class__.__name__, self.name)
