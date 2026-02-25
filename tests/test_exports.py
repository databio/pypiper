"""Verify pypiper's public API exports are complete and correct."""

import pypiper

# Every symbol that downstream pipelines depend on
REQUIRED_PUBLIC_API = [
    "PipelineManager",
    "NGSTk",
    "NGSTools",
    "Pipeline",
    "Stage",
    "add_pypiper_args",
    "build_command",
    "check_all_commands",
    "determine_uncallable",
    "get_first_value",
    "head",
    "is_fastq",
    "is_gzipped_fastq",
    "is_unzipped_fastq",
    "is_sam_or_bam",
    "logger_via_cli",
    "result_formatter_markdown",
    "FLAGS",
    "COMPLETE_FLAG",
    "FAIL_FLAG",
    "RUN_FLAG",
    "WAIT_FLAG",
    "PAUSE_FLAG",
]


class TestPublicAPI:
    def test_required_symbols_accessible(self):
        for name in REQUIRED_PUBLIC_API:
            assert hasattr(pypiper, name), f"pypiper.{name} is not accessible"

    def test_no_private_leaks(self):
        public_names = [n for n in dir(pypiper) if not n.startswith("_")]
        for name in public_names:
            obj = getattr(pypiper, name)
            if callable(obj) and hasattr(obj, "__name__"):
                assert not obj.__name__.startswith("_"), (
                    f"Private function {obj.__name__} leaked into pypiper namespace"
                )
