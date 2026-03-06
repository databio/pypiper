"""Tests for pipestat record_identifier, stats_dict, and pipeline_stats_file fixes.

Covers:
- from_config() path sets record_identifier unconditionally
- report_result() populates stats_dict
- get_stat() retrieves values reported via report_result()
- pipeline_stats_file matches pipestat file in config path
"""

import os

import yaml

from pypiper import PipelineManager

TEST_SCHEMA_PATH = os.path.join(
    os.path.dirname(__file__), "data", "test_pipestat_output_schema.yaml"
)


def _write_pipestat_config(tmpdir, schema_path, results_file_path=None):
    """Write a minimal pipestat config YAML and return its path."""
    if results_file_path is None:
        results_file_path = os.path.join(str(tmpdir), "results.yaml")
    config = {
        "schema_path": schema_path,
        "results_file_path": results_file_path,
    }
    config_path = os.path.join(str(tmpdir), "pipestat_config.yaml")
    with open(config_path, "w") as f:
        yaml.dump(config, f)
    return config_path


class TestFromConfigRecordIdentifier:
    """Test that from_config() path always sets record_identifier from pypiper."""

    def test_record_identifier_set_from_pypiper(self, tmpdir):
        """record_identifier on PipestatManager should match what pypiper provides."""
        config_path = _write_pipestat_config(tmpdir, TEST_SCHEMA_PATH)
        outfolder = os.path.join(str(tmpdir), "output")
        os.makedirs(outfolder, exist_ok=True)

        pm = PipelineManager(
            name="test_pipeline",
            outfolder=outfolder,
            pipestat_config=config_path,
            pipestat_record_identifier="my_sample",
        )
        try:
            assert pm.pipestat.record_identifier == "my_sample"
        finally:
            pm.stop_pipeline()

    def test_record_identifier_overrides_config_value(self, tmpdir):
        """Even if the pipestat config has a record_identifier, pypiper should override it."""
        results_file = os.path.join(str(tmpdir), "results.yaml")
        config = {
            "schema_path": TEST_SCHEMA_PATH,
            "results_file_path": results_file,
            "record_identifier": "stale_record",
        }
        config_path = os.path.join(str(tmpdir), "pipestat_config.yaml")
        with open(config_path, "w") as f:
            yaml.dump(config, f)

        outfolder = os.path.join(str(tmpdir), "output")
        os.makedirs(outfolder, exist_ok=True)

        pm = PipelineManager(
            name="test_pipeline",
            outfolder=outfolder,
            pipestat_config=config_path,
            pipestat_record_identifier="correct_sample",
        )
        try:
            assert pm.pipestat.record_identifier == "correct_sample"
        finally:
            pm.stop_pipeline()

    def test_record_identifier_defaults_to_default_sample_name(self, tmpdir):
        """When no record_identifier is provided, it should default to DEFAULT_SAMPLE_NAME."""
        config_path = _write_pipestat_config(tmpdir, TEST_SCHEMA_PATH)
        outfolder = os.path.join(str(tmpdir), "output")
        os.makedirs(outfolder, exist_ok=True)

        pm = PipelineManager(
            name="test_pipeline",
            outfolder=outfolder,
            pipestat_config=config_path,
        )
        try:
            # DEFAULT_SAMPLE_NAME is "default_pipeline_name"
            assert pm.pipestat.record_identifier is not None
            assert pm.pipestat.record_identifier != ""
        finally:
            pm.stop_pipeline()


class TestReportResultStatsDict:
    """Test that report_result() populates stats_dict and get_stat() works."""

    def test_report_result_populates_stats_dict(self, tmpdir):
        """report_result() should store the value in stats_dict."""
        config_path = _write_pipestat_config(tmpdir, TEST_SCHEMA_PATH)
        outfolder = os.path.join(str(tmpdir), "output")
        os.makedirs(outfolder, exist_ok=True)

        pm = PipelineManager(
            name="test_pipeline",
            outfolder=outfolder,
            pipestat_config=config_path,
            pipestat_record_identifier="sample1",
        )
        try:
            pm.report_result("key1", "value1")
            assert pm.stats_dict["key1"] == "value1"
        finally:
            pm.stop_pipeline()

    def test_get_stat_retrieves_reported_value(self, tmpdir):
        """get_stat() should return a value that was reported via report_result()."""
        config_path = _write_pipestat_config(tmpdir, TEST_SCHEMA_PATH)
        outfolder = os.path.join(str(tmpdir), "output")
        os.makedirs(outfolder, exist_ok=True)

        pm = PipelineManager(
            name="test_pipeline",
            outfolder=outfolder,
            pipestat_config=config_path,
            pipestat_record_identifier="sample1",
        )
        try:
            pm.report_result("key1", "hello")
            result = pm.get_stat("key1")
            assert result == "hello"
        finally:
            pm.stop_pipeline()


class TestPipelineStatsFileSync:
    """Test that pipeline_stats_file is synced to pipestat's file when using config."""

    def test_pipeline_stats_file_matches_pipestat_file(self, tmpdir):
        """When using pipestat_config, pipeline_stats_file should point to pipestat's file."""
        results_file = os.path.join(str(tmpdir), "my_results.yaml")
        config_path = _write_pipestat_config(
            tmpdir, TEST_SCHEMA_PATH, results_file_path=results_file
        )
        outfolder = os.path.join(str(tmpdir), "output")
        os.makedirs(outfolder, exist_ok=True)

        pm = PipelineManager(
            name="test_pipeline",
            outfolder=outfolder,
            pipestat_config=config_path,
            pipestat_record_identifier="sample1",
        )
        try:
            assert pm.pipeline_stats_file == pm.pipestat.file
        finally:
            pm.stop_pipeline()
