"""Tests for _LogTee and logging lifecycle bugs found in PR #239 review.

Each test targets one bug:
1. _LogTee.write() must return character count (TextIOBase contract)
2. Subprocess output must appear exactly once in the log (not doubled)
3. sys.stdout must be restored after a halted pipeline
4. FileHandler must be cleaned up after stop_pipeline()
"""

import logging
import os
import sys

import pypiper
from pypiper.manager import _LogTee

TEST_SCHEMA_PATH = os.path.join(
    os.path.dirname(__file__), "data", "test_pipestat_output_schema.yaml"
)


class TestLogTeeWriteContract:
    """_LogTee.write() must return character count per TextIOBase contract."""

    def test_write_returns_character_count(self, tmp_path):
        """print() through _LogTee works correctly (write returns int)."""
        log_file = str(tmp_path / "test.log")
        original = sys.stdout
        tee = _LogTee(original, log_file)
        result = tee.write("hello")
        assert isinstance(result, int), f"write() returned {type(result).__name__}, expected int"
        assert result == 5

    def test_writelines_mirrors_to_log(self, tmp_path):
        """writelines() copies all lines to the log file."""
        log_file = str(tmp_path / "test.log")
        original = sys.stdout
        tee = _LogTee(original, log_file)
        tee.writelines(["line1\n", "line2\n"])
        with open(log_file) as f:
            content = f.read()
        assert "line1" in content
        assert "line2" in content


class TestSubprocessOutputNotDoubled:
    """Subprocess output must appear exactly once in the pipeline log."""

    def test_subprocess_output_appears_once(self, tmp_path):
        """Command output appears exactly once in the log, not doubled."""
        outfolder = str(tmp_path / "out")
        pm = pypiper.PipelineManager(
            "test_pipe",
            outfolder=outfolder,
            pipestat_schema=TEST_SCHEMA_PATH,
        )
        pm.run("echo UNIQUE_MARKER_ONCE_ONLY", lock_name="t")
        pm.stop_pipeline()

        log_file = os.path.join(outfolder, "test_pipe_log.md")
        with open(log_file) as f:
            lines = f.readlines()
        # Count lines that are the raw output (not the command info line with backticks)
        output_lines = [
            line for line in lines if "UNIQUE_MARKER_ONCE_ONLY" in line and "`" not in line
        ]
        assert len(output_lines) == 1, (
            f"Expected subprocess output once in log, found {len(output_lines)} times"
        )


class TestStdoutRestoredAfterHalt:
    """sys.stdout must be restored even when a pipeline is halted."""

    def test_stdout_restored_after_halt(self, tmp_path):
        """After a halted pipeline, sys.stdout is the original stream."""
        original_stdout = sys.stdout
        outfolder = str(tmp_path / "out")
        pm = pypiper.PipelineManager(
            "test_pipe",
            outfolder=outfolder,
            pipestat_schema=TEST_SCHEMA_PATH,
        )
        pm.halt(raise_error=False)
        # After halt (which calls stop_pipeline), stdout should be restored
        assert sys.stdout is original_stdout, (
            f"sys.stdout is {type(sys.stdout).__name__}, expected original stream"
        )


class TestFileHandlerCleanup:
    """FileHandler must be removed from logger after stop_pipeline()."""

    def test_filehandler_removed_after_stop(self, tmp_path):
        """After stop_pipeline(), the FileHandler is removed from the logger."""
        outfolder = str(tmp_path / "out")
        pm = pypiper.PipelineManager(
            "test_pipe",
            outfolder=outfolder,
            pipestat_schema=TEST_SCHEMA_PATH,
        )
        pm.stop_pipeline()
        file_handlers = [h for h in pm._logger.handlers if isinstance(h, logging.FileHandler)]
        assert len(file_handlers) == 0, (
            f"{len(file_handlers)} FileHandler(s) still on logger after stop"
        )
