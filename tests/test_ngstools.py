"""Tests for NGSTools class."""

from pypiper.ngstk import NGSTk, NGSTools


class TestNGSTools:
    """Test the NGSTools class functionality."""

    def test_default_echo_behavior(self):
        """Test that missing tools return their own name (echo behavior)."""
        tools = NGSTools()
        assert tools.samtools == "samtools"
        assert tools.java == "java"
        assert tools.picard == "picard"
        assert tools.bedtools == "bedtools"

    def test_config_override(self):
        """Test that configured tools use the provided paths."""
        tools = NGSTools({"samtools": "/usr/local/bin/samtools", "java": "/usr/bin/java"})
        assert tools.samtools == "/usr/local/bin/samtools"
        assert tools.java == "/usr/bin/java"
        # Non-configured tools still echo
        assert tools.picard == "picard"
        assert tools.bedtools == "bedtools"

    def test_all_annotations_present(self):
        """Test that all 18 expected tool annotations exist."""
        expected_tools = {
            "samtools",
            "java",
            "picard",
            "bedtools",
            "bowtie2",
            "tophat",
            "sambamba",
            "kallisto",
            "macs2",
            "Rscript",
            "spp",
            "fastqc",
            "skewer",
            "trimmomatic",
            "genomeCoverageBed",
            "bedGraphToBigWig",
            "python",
            "scripts_dir",
        }
        assert set(NGSTools.__annotations__.keys()) == expected_tools
        assert len(NGSTools.__annotations__) == 18

    def test_all_tools_initialized(self):
        """Test that all annotated tools are initialized."""
        tools = NGSTools()
        for attr in NGSTools.__annotations__:
            assert hasattr(tools, attr)
            assert getattr(tools, attr) == attr  # echo behavior


class TestNGSTkWithNGSTools:
    """Test NGSTk integration with NGSTools."""

    def test_ngstk_standalone(self):
        """Test NGSTk with no pm or config_file."""
        tk = NGSTk()
        assert isinstance(tk.tools, NGSTools)
        assert tk.tools.samtools == "samtools"
        assert tk.tools.java == "java"

    def test_ngstk_command_building(self):
        """Test that NGSTk builds commands correctly with default tools."""
        tk = NGSTk()
        cmd = tk.samtools_index("test.bam")
        assert cmd == "samtools index test.bam"

    def test_ngstk_with_custom_tools(self):
        """Test NGSTk with a mock pm that has custom tool paths."""

        class MockConfig:
            def __init__(self):
                self.tools = {"samtools": "/opt/samtools/bin/samtools"}
                self.parameters = {}

        class MockPM:
            def __init__(self):
                self.config = MockConfig()
                self.cores = 1

        pm = MockPM()
        tk = NGSTk(pm=pm)
        assert tk.tools.samtools == "/opt/samtools/bin/samtools"
        # Non-configured tools still echo
        assert tk.tools.java == "java"

    def test_ngstk_parameters_are_dict(self):
        """Test that parameters is a plain dict."""
        tk = NGSTk()
        assert isinstance(tk.parameters, dict)
        assert tk.parameters == {}

    def test_ngstk_with_none_config(self):
        """Test NGSTk handles pm with None config gracefully."""

        class MockPM:
            def __init__(self):
                self.config = None
                self.cores = 1

        pm = MockPM()
        tk = NGSTk(pm=pm)
        assert isinstance(tk.tools, NGSTools)
        assert tk.tools.samtools == "samtools"
