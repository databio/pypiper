"""Characterization tests for logging and subprocess output capture.

These tests pin the required behavior: subprocess stdout/stderr must appear
in the log file, piped commands must capture stderr from all stages,
and multiple managers must produce independent logs.

Tests run pipeline scripts as subprocesses (via temp script files, not -c)
to test real tee-mode pipeline behavior. After the tee-removal refactoring,
these same assertions must still hold.
"""

import os
import subprocess
import sys
import tempfile
import textwrap

TEST_SCHEMA_PATH = os.path.join(
    os.path.dirname(__file__), "data", "test_pipestat_output_schema.yaml"
)


def _run_pipeline_script(script: str, timeout: int = 30) -> subprocess.CompletedProcess:
    """Write script to a temp file and run it, so __main__.__file__ is set."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".py", delete=False) as f:
        f.write(script)
        script_path = f.name
    try:
        result = subprocess.run(
            [sys.executable, script_path],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    finally:
        os.unlink(script_path)
    return result


class TestSubprocessOutputCapture:
    """Subprocess stdout and stderr must appear in the log file."""

    def test_stdout_in_log(self, tmp_path):
        """Command stdout appears in the log file."""
        outfolder = str(tmp_path / "out")
        script = textwrap.dedent(f"""\
            import pypiper
            pm = pypiper.PipelineManager(
                "test_pipe", outfolder="{outfolder}",
                pipestat_schema="{TEST_SCHEMA_PATH}",
            )
            pm.run("echo MARKER_STDOUT_12345", lock_name="t")
            pm.stop_pipeline()
        """)
        result = _run_pipeline_script(script)
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        log_file = os.path.join(outfolder, "test_pipe_log.md")
        assert os.path.isfile(log_file), f"Log file not created: {log_file}"
        with open(log_file) as f:
            log_content = f.read()
        assert "MARKER_STDOUT_12345" in log_content

    def test_stderr_in_log(self, tmp_path):
        """Command stderr appears in the log file."""
        outfolder = str(tmp_path / "out")
        script = textwrap.dedent(f"""\
            import pypiper
            pm = pypiper.PipelineManager(
                "test_pipe", outfolder="{outfolder}",
                pipestat_schema="{TEST_SCHEMA_PATH}",
            )
            pm.run("echo MARKER_STDERR_67890 >&2", lock_name="t", shell=True)
            pm.stop_pipeline()
        """)
        result = _run_pipeline_script(script)
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        log_file = os.path.join(outfolder, "test_pipe_log.md")
        with open(log_file) as f:
            log_content = f.read()
        assert "MARKER_STDERR_67890" in log_content

    def test_stdout_and_stderr_both_in_log(self, tmp_path):
        """Both stdout and stderr from the same command appear in the log."""
        outfolder = str(tmp_path / "out")
        script = textwrap.dedent(f"""\
            import pypiper
            pm = pypiper.PipelineManager(
                "test_pipe", outfolder="{outfolder}",
                pipestat_schema="{TEST_SCHEMA_PATH}",
            )
            pm.run(
                "echo OUT_MARKER_111 && echo ERR_MARKER_222 >&2",
                lock_name="t", shell=True,
            )
            pm.stop_pipeline()
        """)
        result = _run_pipeline_script(script)
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        log_file = os.path.join(outfolder, "test_pipe_log.md")
        with open(log_file) as f:
            log_content = f.read()
        assert "OUT_MARKER_111" in log_content
        assert "ERR_MARKER_222" in log_content


class TestPipedCommandCapture:
    """Piped commands must capture stderr from all stages."""

    def test_piped_stderr_from_first_stage(self, tmp_path):
        """stderr from the first process in a pipe is captured in the log."""
        outfolder = str(tmp_path / "out")
        script = textwrap.dedent(f"""\
            import pypiper
            pm = pypiper.PipelineManager(
                "test_pipe", outfolder="{outfolder}",
                pipestat_schema="{TEST_SCHEMA_PATH}",
            )
            pm.run(
                "bash -c 'echo PIPE_ERR_FIRST >&2; echo data' | cat",
                lock_name="t", shell=True,
            )
            pm.stop_pipeline()
        """)
        result = _run_pipeline_script(script)
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        log_file = os.path.join(outfolder, "test_pipe_log.md")
        with open(log_file) as f:
            log_content = f.read()
        assert "PIPE_ERR_FIRST" in log_content

    def test_piped_stdout_from_last_stage(self, tmp_path):
        """stdout from the last process in a pipe appears in the log."""
        outfolder = str(tmp_path / "out")
        script = textwrap.dedent(f"""\
            import pypiper
            pm = pypiper.PipelineManager(
                "test_pipe", outfolder="{outfolder}",
                pipestat_schema="{TEST_SCHEMA_PATH}",
            )
            pm.run("echo PIPE_DATA | cat", lock_name="t")
            pm.stop_pipeline()
        """)
        result = _run_pipeline_script(script)
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        log_file = os.path.join(outfolder, "test_pipe_log.md")
        with open(log_file) as f:
            log_content = f.read()
        assert "PIPE_DATA" in log_content


class TestMultipleManagers:
    """Multiple PipelineManagers in one process produce independent logs."""

    def test_sequential_managers_independent_logs(self, tmp_path):
        """Two managers produce correct, independent log files."""
        out1 = str(tmp_path / "out1")
        out2 = str(tmp_path / "out2")
        script = textwrap.dedent(f"""\
            import pypiper
            pm1 = pypiper.PipelineManager(
                "mgr_one", outfolder="{out1}",
                pipestat_schema="{TEST_SCHEMA_PATH}",
            )
            pm1.run("echo UNIQUE_PM1_AAA", lock_name="t1")
            pm1.stop_pipeline()

            pm2 = pypiper.PipelineManager(
                "mgr_two", outfolder="{out2}",
                pipestat_schema="{TEST_SCHEMA_PATH}",
            )
            pm2.run("echo UNIQUE_PM2_BBB", lock_name="t2")
            pm2.stop_pipeline()
        """)
        result = _run_pipeline_script(script)
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        log1 = os.path.join(out1, "mgr_one_log.md")
        log2 = os.path.join(out2, "mgr_two_log.md")

        assert os.path.isfile(log1), f"Log file 1 not created: {log1}"
        assert os.path.isfile(log2), f"Log file 2 not created: {log2}"

        with open(log1) as f:
            content1 = f.read()
        with open(log2) as f:
            content2 = f.read()

        assert "UNIQUE_PM1_AAA" in content1
        assert "UNIQUE_PM2_BBB" in content2
        # Each log should NOT contain the other's marker
        assert "UNIQUE_PM2_BBB" not in content1
        assert "UNIQUE_PM1_AAA" not in content2


class TestPypiperMessagesInLog:
    """Pypiper's own messages (info, timestamp) must appear in the log file."""

    def test_timestamp_in_log(self, tmp_path):
        """pm.timestamp() messages appear in the log file."""
        outfolder = str(tmp_path / "out")
        script = textwrap.dedent(f"""\
            import pypiper
            pm = pypiper.PipelineManager(
                "test_pipe", outfolder="{outfolder}",
                pipestat_schema="{TEST_SCHEMA_PATH}",
            )
            pm.timestamp("TIMESTAMP_MARKER_XYZ")
            pm.stop_pipeline()
        """)
        result = _run_pipeline_script(script)
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        log_file = os.path.join(outfolder, "test_pipe_log.md")
        with open(log_file) as f:
            log_content = f.read()
        assert "TIMESTAMP_MARKER_XYZ" in log_content

    def test_info_in_log(self, tmp_path):
        """pm.info() messages appear in the log file."""
        outfolder = str(tmp_path / "out")
        script = textwrap.dedent(f"""\
            import pypiper
            pm = pypiper.PipelineManager(
                "test_pipe", outfolder="{outfolder}",
                pipestat_schema="{TEST_SCHEMA_PATH}",
            )
            pm.info("INFO_MARKER_789")
            pm.stop_pipeline()
        """)
        result = _run_pipeline_script(script)
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        log_file = os.path.join(outfolder, "test_pipe_log.md")
        with open(log_file) as f:
            log_content = f.read()
        assert "INFO_MARKER_789" in log_content


class TestPrintCapture:
    """Bare print() calls must appear in the log file."""

    def test_print_captured_in_log(self, tmp_path):
        """Bare print() calls between pm.run() appear in the log file."""
        outfolder = str(tmp_path / "out")
        script = textwrap.dedent(f"""\
            import pypiper
            pm = pypiper.PipelineManager(
                "test_pipe", outfolder="{outfolder}",
                pipestat_schema="{TEST_SCHEMA_PATH}",
            )
            print("PRINT_MARKER_ABC123")
            pm.stop_pipeline()
        """)
        result = _run_pipeline_script(script)
        assert result.returncode == 0, f"Script failed: {result.stderr}"
        log_file = os.path.join(outfolder, "test_pipe_log.md")
        with open(log_file) as f:
            log_content = f.read()
        assert "PRINT_MARKER_ABC123" in log_content


class TestMultiModeLogging:
    """Tests for logging via FileHandler and per-command capture."""

    def test_multi_mode_stdout_in_log(self, tmp_path):
        """In multi mode, subprocess stdout should still appear in the log."""
        outfolder = str(tmp_path / "out")
        script = textwrap.dedent(f"""\
            import pypiper
            pm = pypiper.PipelineManager(
                "test_pipe", outfolder="{outfolder}",
                pipestat_schema="{TEST_SCHEMA_PATH}",
            )
            pm.run("echo MULTI_STDOUT_MARKER", lock_name="t")
            pm.stop_pipeline()
        """)
        result = _run_pipeline_script(script)
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        log_file = os.path.join(outfolder, "test_pipe_log.md")
        assert os.path.isfile(log_file), f"Log file not created: {log_file}"
        with open(log_file) as f:
            log_content = f.read()
        assert "MULTI_STDOUT_MARKER" in log_content

    def test_multi_mode_info_in_log(self, tmp_path):
        """In multi mode, pypiper messages should still appear in the log."""
        outfolder = str(tmp_path / "out")
        script = textwrap.dedent(f"""\
            import pypiper
            pm = pypiper.PipelineManager(
                "test_pipe", outfolder="{outfolder}",
                pipestat_schema="{TEST_SCHEMA_PATH}",
            )
            pm.info("MULTI_INFO_MARKER")
            pm.stop_pipeline()
        """)
        result = _run_pipeline_script(script)
        assert result.returncode == 0, f"Script failed: {result.stderr}"

        log_file = os.path.join(outfolder, "test_pipe_log.md")
        assert os.path.isfile(log_file), f"Log file not created: {log_file}"
        with open(log_file) as f:
            log_content = f.read()
        assert "MULTI_INFO_MARKER" in log_content
