"""Tests for actionable error messages.

These are string-content assertions that verify key error messages contain
actionable guidance. They prevent regressions where someone shortens an
error message back to a non-actionable form.
"""

from pypiper.exceptions import (
    MissingCheckpointError,
    SubprocessError,
    UnknownPipelineStageError,
)


class TestSubprocessError:
    """SubprocessError stores returncode and cmd as structured data."""

    def test_stores_returncode(self):
        e = SubprocessError("failed", returncode=1)
        assert e.returncode == 1

    def test_stores_cmd(self):
        e = SubprocessError("failed", cmd="samtools sort")
        assert e.cmd == "samtools sort"

    def test_defaults_to_none(self):
        e = SubprocessError("failed")
        assert e.returncode is None
        assert e.cmd is None

    def test_message_preserved(self):
        e = SubprocessError("something went wrong")
        assert str(e) == "something went wrong"


class TestMissingCheckpointError:
    """MissingCheckpointError message includes actionable hints."""

    def test_contains_new_start_hint(self):
        e = MissingCheckpointError("align", "/out/align.checkpoint")
        msg = str(e)
        assert "new_start" in msg

    def test_contains_checkpoint_name(self):
        e = MissingCheckpointError("align", "/out/align.checkpoint")
        msg = str(e)
        assert "align" in msg

    def test_contains_filepath(self):
        e = MissingCheckpointError("align", "/out/align.checkpoint")
        msg = str(e)
        assert "/out/align.checkpoint" in msg

    def test_contains_start_point_hint(self):
        e = MissingCheckpointError("align", "/out/align.checkpoint")
        msg = str(e)
        assert "start_point" in msg


class TestUnknownPipelineStageError:
    """UnknownPipelineStageError lists available stages when pipeline is provided."""

    def test_contains_stage_name(self):
        e = UnknownPipelineStageError("nonexistent")
        msg = str(e)
        assert "nonexistent" in msg

    def test_contains_unknown_label(self):
        e = UnknownPipelineStageError("nonexistent")
        msg = str(e)
        assert "Unknown pipeline stage" in msg

    def test_lists_available_stages_with_pipeline(self):
        class FakePipeline:
            def stages(self):
                return ["trim", "align", "count"]

        e = UnknownPipelineStageError("nonexistent", FakePipeline())
        msg = str(e)
        assert "trim" in msg
        assert "align" in msg
        assert "count" in msg
        assert "Available stages" in msg

    def test_contains_typo_hint_with_pipeline(self):
        class FakePipeline:
            def stages(self):
                return ["trim", "align"]

        e = UnknownPipelineStageError("trm", FakePipeline())
        msg = str(e)
        assert "typos" in msg.lower() or "Check for typos" in msg
