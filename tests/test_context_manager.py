"""Tests for PipelineManager context manager support."""

import os

import pytest

from pypiper import PipelineManager


class TestContextManager:
    """Tests for using PipelineManager as a context manager."""

    def test_context_manager_clean_exit(self, tmpdir):
        """On clean exit, stop_pipeline() is called and status is completed."""
        with PipelineManager("test_cm", tmpdir.strpath) as pm:
            pass
        assert pm._has_exit_status

    def test_context_manager_exception_sets_failed(self, tmpdir):
        """Exception inside with block calls fail_pipeline() and propagates."""
        with pytest.raises(ValueError, match="test error"):
            with PipelineManager("test_cm_fail", tmpdir.strpath) as pm:
                raise ValueError("test error")
        assert pm._failed

    def test_context_manager_yields_self(self, tmpdir):
        """The yielded object is the PipelineManager instance."""
        with PipelineManager("test_cm_self", tmpdir.strpath) as pm:
            assert isinstance(pm, PipelineManager)
            pm_ref = pm
        assert pm_ref is pm

    def test_context_manager_run_command(self, tmpdir):
        """End-to-end: run a command with target inside with, verify target created."""
        target = os.path.join(tmpdir.strpath, "output.txt")
        with PipelineManager("test_cm_run", tmpdir.strpath) as pm:
            pm.run(f"echo 'hello' > {target}", target=target)
        assert os.path.isfile(target)

    def test_context_manager_keyboard_interrupt(self, tmpdir):
        """KeyboardInterrupt triggers fail_pipeline() and propagates."""
        with pytest.raises(KeyboardInterrupt):
            with PipelineManager("test_cm_kb", tmpdir.strpath) as pm:
                raise KeyboardInterrupt()
        assert pm._failed
