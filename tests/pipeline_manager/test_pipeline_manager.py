#!/usr/bin/env python
"""PipelineManager tests."""

import glob
import os
import shutil
import subprocess
import time

import pytest

import pypiper
from pypiper.exceptions import SubprocessError
from pypiper.utils import pipeline_filepath

__author__ = "Nathan Sheffield"
__email__ = "nathan@code.databio.org"

# Path to the test schema for pipestat
TEST_SCHEMA_PATH = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    "Data",
    "test_pipestat_output_schema.yaml"
)


@pytest.fixture
def single_pipeline_manager(tmpdir, request):
    """
    Create a single pipeline manager for tests that only need one.

    More efficient than pipeline_managers for tests that don't need
    multiple managers sharing output folders.
    """
    outfolder = str(tmpdir.mkdir("pipeline_output"))
    pp = pypiper.PipelineManager(
        "sample_pipeline", outfolder=outfolder, multi=True,
        pipestat_schema=TEST_SCHEMA_PATH
    )

    def cleanup():
        """Clean up pipeline and clear status flags."""
        try:
            pp.pipestat.clear_status("sample_pipeline", flag_names=["failed"])
        except Exception:
            pass
        try:
            if not pp.status:
                pp.stop_pipeline()
        except Exception:
            pass

    request.addfinalizer(cleanup)
    return pp


@pytest.fixture
def pipeline_managers(tmpdir, request):
    """
    Create three pipeline managers for testing.

    Returns a tuple of (pp, pp2, pp3) where pp and pp2 share the same
    output folder and pp3 has a separate output folder.

    Use single_pipeline_manager for tests that only need one manager.
    """
    outfolder = str(tmpdir.mkdir("pipeline_output"))
    outfolder3 = str(tmpdir.mkdir("pipeline_output3"))

    pp = pypiper.PipelineManager(
        "sample_pipeline", outfolder=outfolder, multi=True,
        pipestat_schema=TEST_SCHEMA_PATH
    )
    pp2 = pypiper.PipelineManager(
        "sample_pipeline2", outfolder=outfolder, multi=True,
        pipestat_schema=TEST_SCHEMA_PATH
    )
    pp3 = pypiper.PipelineManager(
        "sample_pipeline3", outfolder=outfolder3, multi=True,
        pipestat_schema=TEST_SCHEMA_PATH
    )

    def cleanup():
        """Clean up pipelines and clear status flags."""
        # Clear status flags first (before stopping pipelines)
        for pm, name in [(pp, "sample_pipeline"), (pp2, "sample_pipeline2"), (pp3, "sample_pipeline3")]:
            try:
                pm.pipestat.clear_status(name, flag_names=["failed"])
            except Exception:
                pass  # Ignore errors during cleanup

        # Stop all pipelines - only if not already stopped
        for pm in [pp, pp2, pp3]:
            try:
                # Check if pipeline has already been stopped via the 'status' property
                if not pm.status:
                    pm.stop_pipeline()
            except Exception:
                pass  # Ignore errors during cleanup

    request.addfinalizer(cleanup)
    return pp, pp2, pp3


def _is_file(pp, filename):
    """Determine if the pipeline manager has this file."""
    filepath = pipeline_filepath(pp, filename=filename)
    return os.path.isfile(filepath)


def _assert_file(pp, filename):
    """Assert that the named file exists for pipeline manager."""
    if not _is_file(pp, filename):
        outfolder_contents = os.listdir(pp.outfolder)
        pytest.fail(
            f"Expected file {filename} not found. "
            f"Pipeline outfolder contents:\n{chr(10).join(outfolder_contents)}"
        )


def _assert_not_file(pp, filename):
    """Assert that given file doesn't exist for manager."""
    assert not _is_file(pp, filename)


def _assert_lines(expected, observed):
    """Assert equality between collections of lines."""
    if isinstance(observed, str) and os.path.isfile(observed):
        with open(observed, "r") as f:
            observed = f.readlines()
    assert expected == [line.strip() for line in observed]


class TestPipelineManagerInitialization:
    """Tests for PipelineManager initialization."""

    def test_pipeline_names(self, pipeline_managers):
        """Test that pipeline managers have correct names."""
        pp, pp2, pp3 = pipeline_managers
        assert pp.name == "sample_pipeline"
        assert pp2.name == "sample_pipeline2"
        assert pp3.name == "sample_pipeline3"

    def test_outfolder_creation(self, pipeline_managers):
        """Test that output folder is created."""
        pp, pp2, pp3 = pipeline_managers
        assert os.path.isdir(pp.outfolder)
        assert os.path.isdir(pp3.outfolder)


class TestStatusFlags:
    """Tests for status flag functionality."""

    def test_status_flag_completed(self, single_pipeline_manager):
        """Test setting completed status flag."""
        pp = single_pipeline_manager
        pp._set_status_flag("completed")
        _assert_file(pp, "sample_pipeline_DEFAULT_SAMPLE_NAME_completed.flag")

    def test_status_flag_running(self, single_pipeline_manager):
        """Test setting running status flag clears completed."""
        pp = single_pipeline_manager
        pp._set_status_flag("completed")
        pp._set_status_flag("running")
        _assert_not_file(pp, "sample_pipeline_DEFAULT_SAMPLE_NAME_testing.flag")
        _assert_file(pp, "sample_pipeline_DEFAULT_SAMPLE_NAME_running.flag")


class TestLockWaiting:
    """Tests for lock waiting functionality."""

    def test_wait_for_lock(self, pipeline_managers):
        """Test that pipeline waits for locks."""
        pp, pp2, pp3 = pipeline_managers
        pp2.wait = False
        pp.wait = False
        sleep_lock = pipeline_filepath(pp, filename="lock.sleep")
        # Use short sleep (0.1s) - just enough to prove lock waiting works
        subprocess.Popen("sleep .1; rm " + sleep_lock, shell=True)
        pp._create_file(sleep_lock)
        cmd = "echo hello"
        stamp = time.time()
        pp.run(cmd, lock_name="sleep")
        # Verify pipeline waited (at least one 0.5s sleep iteration)
        assert pp.time_elapsed(stamp) > 0.4

        # Wait for subprocess
        for p in pp.running_procs.copy():
            pp._wait_for_process(pp.running_procs[p]["p"])
        pp2.wait = True
        pp.wait = True


class TestTargetRespecting:
    """Tests for respecting existing target files."""

    def test_respects_existing_files(self, single_pipeline_manager):
        """Test that pipeline respects files already existing."""
        pp = single_pipeline_manager
        target = pipeline_filepath(pp, filename="tgt")
        if os.path.isfile(target):
            os.remove(target)

        pp.run("echo first > " + target, target, shell=True)
        # Should not run
        pp.run("echo second > " + target, target, shell=True)
        with open(target) as f:
            lines = f.readlines()
        _assert_lines(["first"], lines)

    def test_targetless_command(self, single_pipeline_manager):
        """Test executing a targetless command."""
        pp = single_pipeline_manager
        target = pipeline_filepath(pp, filename="tgt")
        pp.run("echo third > " + target, target=None, lock_name="test", shell=True)
        with open(target) as f:
            lines = f.readlines()
        _assert_lines(["third"], lines)


class TestResultReporting:
    """Tests for result reporting functionality."""

    def test_report_and_get_result(self, single_pipeline_manager):
        """Test reporting and retrieving results."""
        pp = single_pipeline_manager
        pp.report_result("key1", "abc")
        pp.report_result("key2", "def")
        key1 = pp.get_stat("key1")
        assert key1 == "abc"

    def test_result_not_shared_between_pipelines(self, pipeline_managers):
        """Test that results are not shared between different pipelines."""
        pp, pp2, pp3 = pipeline_managers
        pp.report_result("key1", "abc")
        try:
            key1 = pp2.get_stat("key1")
        except KeyError:
            key1 = None
        assert key1 is None


class TestCleanup:
    """Tests for intermediate file cleanup functionality."""

    def test_dirty_mode_no_cleanup(self, single_pipeline_manager):
        """Test that dirty mode prevents file cleanup."""
        pp = single_pipeline_manager
        tgt1 = pipeline_filepath(pp, filename="tgt1.temp")
        tgt2 = pipeline_filepath(pp, filename="tgt2.temp")
        tgt3 = pipeline_filepath(pp, filename="tgt3.temp")
        tgt4 = pipeline_filepath(pp, filename="tgt4.txt")

        pp.run("touch " + tgt1 + " " + tgt2 + " " + tgt3 + " " + tgt4, lock_name="test")

        # In global dirty mode, even non-manual clean files should not be deleted
        pp.dirty = True
        pp.clean_add(tgt3)
        pp.clean_add(pipeline_filepath(pp, filename="*.temp"))
        pp.clean_add(tgt4)
        pp._cleanup()

        assert os.path.isfile(tgt1)
        assert os.path.isfile(tgt2)
        assert os.path.isfile(tgt3)
        assert os.path.isfile(tgt4)

    def test_regular_mode_cleanup(self, single_pipeline_manager):
        """Test that regular mode cleans up files."""
        pp = single_pipeline_manager
        tgt1 = pipeline_filepath(pp, filename="tgt1.temp")
        tgt2 = pipeline_filepath(pp, filename="tgt2.temp")
        tgt3 = pipeline_filepath(pp, filename="tgt3.temp")
        tgt4 = pipeline_filepath(pp, filename="tgt4.txt")

        pp.run("touch " + tgt1 + " " + tgt2 + " " + tgt3 + " " + tgt4, lock_name="test")

        pp.dirty = False
        pp.clean_add(pipeline_filepath(pp, filename="*.temp"))
        pp.clean_add(tgt4)
        pp._cleanup()

        assert not os.path.isfile(tgt1)
        assert not os.path.isfile(tgt2)
        assert not os.path.isfile(tgt3)
        assert not os.path.isfile(tgt4)

    def test_conditional_cleanup(self, pipeline_managers):
        """Test conditional cleanup based on other pipeline status."""
        pp, pp2, pp3 = pipeline_managers
        tgt5 = pipeline_filepath(pp, filename="tgt5.txt")
        tgt8 = pipeline_filepath(pp, filename="tgt8.cond")
        tgt9 = pipeline_filepath(pp, filename="tgt9.cond")

        pp.run("touch " + tgt5 + " " + tgt8 + " " + tgt9, lock_name="test")

        pp.dirty = False
        pp.clean_add(tgt5, conditional=True)
        pp.clean_add(pipeline_filepath(pp, filename="*.cond"), conditional=True)
        pp._cleanup()

        # Conditional delete should not delete while pp2 is running
        assert os.path.isfile(tgt5)
        assert os.path.isfile(tgt8)
        assert os.path.isfile(tgt9)

        # Stopping pp2 should cause conditional files to be deleted
        pp2.stop_pipeline()
        pp._cleanup()
        assert not os.path.isfile(tgt5)
        assert not os.path.isfile(tgt8)
        assert not os.path.isfile(tgt9)

    def test_manual_cleanup_not_auto_deleted(self, single_pipeline_manager):
        """Test that manual cleanup files are not auto-deleted."""
        pp = single_pipeline_manager
        tgt7 = pipeline_filepath(pp, filename="tgt7.txt")
        pp.run("touch " + tgt7, tgt7)
        pp.clean_add(tgt7, manual=True)
        pp.stop_pipeline()
        # Manual clean should not clean even after pipeline stops
        assert os.path.isfile(tgt7)

    def test_auto_cleanup_on_run(self, pipeline_managers):
        """Test automatic cleanup when clean=True is passed to run."""
        pp, pp2, pp3 = pipeline_managers
        tgt10 = pipeline_filepath(pp, filename="tgt10.txt")
        pp.run("touch " + tgt10, target=tgt10, clean=True)
        assert os.path.isfile(tgt10)  # File exists initially

        pp2.stop_pipeline()  # Stop pp2 so conditional cleanup can proceed
        pp._cleanup()
        assert not os.path.isfile(tgt10)

    def test_cleanup_script_contents(self, single_pipeline_manager):
        """Test that cleanup script is populated correctly."""
        pp = single_pipeline_manager
        tgt3 = pipeline_filepath(pp, filename="tgt3.temp")
        tgt6 = pipeline_filepath(pp, filename="tgt6.txt")
        tgt6_abs = os.path.abspath(tgt6)

        pp.run("touch " + tgt3 + " " + tgt6, lock_name="test")

        pp.dirty = True
        pp.clean_add(tgt3)

        # Test from different wd
        cfile = pp.cleanup_file
        ofolder = pp.outfolder
        cwd = os.getcwd()
        pp.clean_add(tgt6_abs)

        os.chdir(pp.outfolder)
        pp.outfolder = "../" + ofolder
        pp.cleanup_file = "../" + cfile
        pp.clean_add(tgt6_abs)
        os.chdir(cwd)
        pp.cleanup_file = cfile
        pp.outfolder = ofolder

        with open(pp.cleanup_file) as f:
            lines = f.readlines()

        assert lines[2] == "rm tgt3.temp\n"

    def test_clean_add_none_is_noop(self, single_pipeline_manager):
        """Test that clean_add(None) silently does nothing."""
        pp = single_pipeline_manager
        pp.clean_add(None)
        pp.clean_add(None, conditional=True)
        pp.clean_add(None, manual=True)
        assert len(pp.cleanup_list) == 0
        assert len(pp.cleanup_list_conditional) == 0

    def test_clean_add_none_among_real_files(self, single_pipeline_manager):
        """Test that None values mixed with real paths work correctly."""
        pp = single_pipeline_manager
        tgt = pipeline_filepath(pp, filename="real_file.txt")
        pp.run("touch " + tgt, lock_name="test")
        pp.clean_add(None)
        pp.clean_add(tgt)
        pp.clean_add(None, conditional=True)
        assert len(pp.cleanup_list) == 1
        assert pp.cleanup_list[0] == tgt
        assert len(pp.cleanup_list_conditional) == 0

    def test_clean_add_none_no_cleanup_file_created(self, single_pipeline_manager):
        """Test that clean_add(None, manual=True) does not create a cleanup script."""
        pp = single_pipeline_manager
        pp.clean_add(None, manual=True)
        assert not os.path.exists(pp.cleanup_file)


class TestFailureHandling:
    """Tests for failure and nofail options."""

    def test_nofail_option(self, single_pipeline_manager):
        """Test that nofail option prevents exception raising."""
        pp = single_pipeline_manager
        cmd = "thiscommandisbad"
        # Should not raise an error
        pp.run(cmd, target=None, lock_name="badcommand", nofail=True)
        pp.callprint(cmd, shell=None, nofail=True)

    def test_failure_raises_exception(self, single_pipeline_manager):
        """Test that bad command raises SubprocessError."""
        pp = single_pipeline_manager
        cmd = "thiscommandisbad"
        with pytest.raises(SubprocessError):
            pp.run(cmd, target=None, lock_name="badcommand")


class TestDynamicRecovery:
    """Tests for dynamic recovery functionality."""

    def test_signal_int_handler(self, single_pipeline_manager):
        """Test signal interrupt handler."""
        pp = single_pipeline_manager
        pp.locks.append("lock.sleep")
        with pytest.raises(KeyboardInterrupt):
            pp._signal_int_handler(None, None)

    def test_recover_from_existing_lock(self, single_pipeline_manager):
        """Test recovery when lock file already exists (overwrite locks mode)."""
        pp = single_pipeline_manager
        # Enable recover mode to overwrite existing locks
        pp.overwrite_locks = True
        sleep_lock = pipeline_filepath(pp, filename="lock.sleep")
        pp._create_file(sleep_lock)
        cmd = "echo hello"
        # This should succeed because overwrite_locks is True
        pp.run(cmd, lock_name="sleep")


class TestNewStart:
    """Tests for new_start functionality."""

    def test_new_start_reruns_commands(self, single_pipeline_manager):
        """Test that new_start causes commands to rerun."""
        pp = single_pipeline_manager
        target = pipeline_filepath(pp, filename="tgt")
        if os.path.isfile(target):
            os.remove(target)

        pp.run("echo first > " + target, target, shell=True)
        pp.run("echo second > " + target, target, shell=True)
        with open(target) as f:
            lines = f.readlines()
        _assert_lines(["first"], lines)

        pp.new_start = True
        # Should run because new_start is True
        pp.run("echo third > " + target, target, shell=True)
        with open(target) as f:
            lines = f.readlines()
        _assert_lines(["third"], lines)


class TestDualTarget:
    """Tests for dual target functionality."""

    def test_single_target_exists_skips(self, single_pipeline_manager):
        """Test that command is skipped if single target exists."""
        pp = single_pipeline_manager
        tgt1 = pipeline_filepath(pp, filename="tgt1.txt")
        tgt6 = pipeline_filepath(pp, filename="tgt6.txt")

        pp.new_start = False
        pp.run("touch " + tgt6, tgt6)
        assert os.path.isfile(tgt6)

        # if target exists, should not run
        pp.run("touch " + tgt1, tgt6)
        assert not os.path.isfile(tgt1)

    def test_dual_target_one_missing_runs(self, single_pipeline_manager):
        """Test that command runs if one of two targets is missing."""
        pp = single_pipeline_manager
        tgt1 = pipeline_filepath(pp, filename="tgt1.txt")
        tgt5 = pipeline_filepath(pp, filename="tgt5.txt")
        tgt6 = pipeline_filepath(pp, filename="tgt6.txt")

        pp.new_start = False
        pp.run("touch " + tgt6, tgt6)

        # if two targets, only one exists, should run
        assert not os.path.isfile(tgt5)
        pp.run("touch " + tgt1, [tgt5, tgt6])
        assert os.path.isfile(tgt1)

    def test_dual_target_both_exist_skips(self, single_pipeline_manager):
        """Test that command is skipped if both targets exist."""
        pp = single_pipeline_manager
        tgt1 = pipeline_filepath(pp, filename="tgt1.txt")
        tgt5 = pipeline_filepath(pp, filename="tgt5.txt")
        tgt6 = pipeline_filepath(pp, filename="tgt6.txt")

        pp.new_start = False
        pp.run("touch " + tgt1 + " " + tgt6, lock_name="setup")

        # if two targets, both exist, should not run
        assert not os.path.isfile(tgt5)
        pp.run("touch " + tgt5, [tgt1, tgt6])
        assert not os.path.isfile(tgt5)


def _make_pipe_filepath(pm, filename):
    return os.path.join(pm.outfolder, filename)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
