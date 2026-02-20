"""Context manager for temporarily changing folder."""

import os
from types import TracebackType

__all__ = []  # Internal, not part of public API


class FolderContext(object):
    """Context manager for temporarily changing directory."""

    def __init__(self, folder: str) -> None:
        """Store the previous working path to restore upon exit.

        Args:
            folder: Path to set as new working directory.
        """
        if not os.path.isdir(folder):
            raise ValueError("Requested temp entry to non-folder: {}".format(folder))
        self._prevdir = os.getcwd()
        self._currdir = folder

    def __enter__(self) -> None:
        """Make the working directory switch."""
        os.chdir(self._currdir)

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        """Switch back to the previous working directory."""
        if not os.path.isdir(self._prevdir):
            raise RuntimeError("Return path is no longer a directory: {}".format(self._prevdir))
        os.chdir(self._prevdir)
