"""Tests verifying that stdlib replacements for pandas produce correct results."""

import textwrap

import pytest

from ubiquerg import parse_timedelta
from pypiper.ngstk import NGSTk


class TestGetElapsedTime:
    """Tests for PipelineManager.get_elapsed_time() using profile TSV fixtures."""

    def test_dedup_keeps_last(self, tmp_path):
        """Profile with duplicate cid keeps last occurrence, sums correctly."""
        profile = tmp_path / "test_profile.tsv"
        profile.write_text(textwrap.dedent("""\
            # Pipeline started at 01-15 10:30:00

            # pid\thash\tcid\truntime\tmem\tcmd\tlock
            12345\tabc123\t1\t0:00:05\t 128.5000\techo hello\tlock.echo_hello
            12345\tabc123\t2\t0:01:30\t 256.7500\tsamtools sort\tlock.samtools_sort
            12345\tabc123\t3\t0:00:45.50\t 64.2500\tfastqc input\tlock.fastqc
            12345\tabc123\t2\t0:02:00\t 512.0000\tsamtools sort\tlock.samtools_sort
        """))

        class FakeManager:
            pipeline_profile_file = str(profile)
            starttime = 0
            def time_elapsed(self, st): return 0.0

        from pypiper.manager import PipelineManager
        mgr = FakeManager()
        result = PipelineManager.get_elapsed_time(mgr)
        # cid 1: 5s, cid 2: 120s (last), cid 3: 45.5s => 170.5
        assert result == pytest.approx(170.5)

    def test_missing_profile_falls_back(self, tmp_path):
        """When profile file does not exist, falls back to time_elapsed()."""
        class FakeManager:
            pipeline_profile_file = str(tmp_path / "nonexistent.tsv")
            starttime = 0
            def time_elapsed(self, st): return 42.0

        from pypiper.manager import PipelineManager
        mgr = FakeManager()
        assert PipelineManager.get_elapsed_time(mgr) == 42.0

    def test_empty_profile_returns_zero(self, tmp_path):
        """Profile with only comments returns 0."""
        profile = tmp_path / "test_profile.tsv"
        profile.write_text("# Pipeline started at 01-01 00:00:00\n\n# pid\thash\tcid\truntime\tmem\tcmd\tlock\n")

        class FakeManager:
            pipeline_profile_file = str(profile)
            starttime = 0
            def time_elapsed(self, st): return 99.0

        from pypiper.manager import PipelineManager
        mgr = FakeManager()
        assert PipelineManager.get_elapsed_time(mgr) == 0


class TestNoPandasImport:
    """Verify that replaced functions do not import pandas."""

    def test_parse_timedelta_no_pandas(self):
        import sys
        parse_timedelta("0:01:00")

    def test_ngstk_methods_no_pandas(self):
        import sys
        was_loaded = "pandas" in sys.modules
        tk = NGSTk()
        tk.parse_bowtie_stats("/nonexistent")
        tk.parse_duplicate_stats("/nonexistent")
        tk.parse_qc("/nonexistent")
        if not was_loaded:
            assert "pandas" not in sys.modules
