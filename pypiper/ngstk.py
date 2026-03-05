"""Optional genomics toolkit with NGS processing convenience functions."""

import errno
import os
import re
import subprocess
from collections.abc import Callable, Iterable
from typing import Any

from yacman import load_yaml

from .exceptions import UnsupportedFiletypeException
from .utils import is_fastq, is_gzipped_fastq, is_sam_or_bam

__all__ = ["NGSTk", "NGSTools"]


class NGSTools:
    """Container for NGS tool paths/commands used by NGSTk.

    Each attribute represents a command-line tool. If a tool path is not
    provided in the config, it defaults to the tool's own name (i.e., assumes
    the tool is on $PATH).

    Example:
        tools = NGSTools({"samtools": "/usr/local/bin/samtools"})
        tools.samtools  # "/usr/local/bin/samtools"
        tools.java      # "java" (not configured, echoes name)
    """

    samtools: str
    java: str
    picard: str
    bedtools: str
    bowtie2: str
    tophat: str
    sambamba: str
    kallisto: str
    macs2: str
    Rscript: str
    spp: str
    fastqc: str
    skewer: str
    trimmomatic: str
    genomeCoverageBed: str
    bedGraphToBigWig: str
    python: str
    scripts_dir: str

    def __init__(self, config: dict[str, str] | None = None) -> None:
        """Initialize tools from config dict.

        Args:
            config: Dict mapping tool names to paths. Missing tools default
                to their own name (echo behavior).
        """
        config = config or {}
        for attr in type(self).__annotations__:
            setattr(self, attr, config.get(attr, attr))


class NGSTk:
    """Build shell command strings for common NGS processing operations.

    Example:
        tk = NGSTk()
        tk.samtools_index("sample.bam")  # => "samtools index sample.bam"

        tk = NGSTk("pipeline_config.yaml")
        tk.samtools_index("sample.bam")  # uses configured samtools path

    Tool paths come from the config's "tools" section; unconfigured tools
    default to their name (assuming they're on $PATH).

    Args:
        config_file: Path to pipeline YAML config file.
        pm: PipelineManager to associate with this toolkit.
    """

    def __init__(self, config_file: str | None = None, pm: Any = None) -> None:
        self.pm = pm

        # Determine tools config from pm, config_file, or empty
        if pm is not None and pm.config is not None:
            # pm.config may be EchoDict or dict - handle both
            if hasattr(pm.config, "tools"):
                tools_cfg = pm.config.tools
                tools_config = dict(tools_cfg) if tools_cfg else {}
            elif isinstance(pm.config, dict):
                tools_config = pm.config.get("tools", {})
            else:
                tools_config = {}
        elif config_file is not None:
            loaded = load_yaml(config_file)
            tools_config = loaded.get("tools", {}) if loaded else {}
        else:
            tools_config = {}

        self.tools = NGSTools(tools_config)

        # Load parameters as a plain dict (no echo needed)
        if pm is not None and pm.config is not None:
            if hasattr(pm.config, "parameters"):
                params = pm.config.parameters
                self.parameters = dict(params) if params else {}
            elif isinstance(pm.config, dict):
                self.parameters = pm.config.get("parameters", {})
            else:
                self.parameters = {}
        else:
            self.parameters = {}

        # If pigz is available, use that. Otherwise, default to gzip.
        if (
            self.pm is not None
            and hasattr(self.pm, "cores")
            and self.pm.cores > 1
            and self.check_command("pigz")
        ):
            self.ziptool_cmd = "pigz -f -p {}".format(self.pm.cores)
        else:
            self.ziptool_cmd = "gzip -f"

    def _ensure_folders(self, *paths: str) -> None:
        """Create parent directories for the given file paths if needed."""
        for p in paths:
            # Only provide assurance for absolute paths.
            if not p or not os.path.isabs(p):
                continue
            # See if what we're assuring is file- or folder-like.
            fpath, fname = os.path.split(p)
            base, ext = os.path.splitext(fname)
            # If there's no extension, ensure that we have the whole path.
            # Otherwise, just ensure that we have path to file's folder.
            self.make_dir(fpath if ext else p)

    @property
    def ziptool(self) -> str:
        """Compression command: 'pigz' if available with multiple cores, else 'gzip'."""
        return self.ziptool_cmd

    def make_dir(self, path: str) -> None:
        """Create directory and all intermediates, no error if exists."""
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    def make_sure_path_exists(self, path: str) -> None:
        """Alias for make_dir"""
        self.make_dir(path)

    # Borrowed from looper
    def check_command(self, command: str) -> bool:
        """Check if a command is callable on the system."""

        # Use `command` to see if command is callable, store exit code
        code = os.system("command -v {0} >/dev/null 2>&1 || {{ exit 1; }}".format(command))

        # If exit code is not 0, report which command failed and return False, else return True
        if code != 0:
            print("Command is not callable: {0}".format(command))
            return False
        else:
            return True

    def get_file_size(self, filenames: str | list[str]) -> float:
        """Get total size of file(s) in megabytes."""
        # use (1024 ** 3) for gigabytes
        # equivalent to: stat -Lc '%s' filename

        # If given a list, recurse through it.
        if type(filenames) is list:
            return sum([self.get_file_size(filename) for filename in filenames])

        return round(
            sum([float(os.stat(f).st_size) for f in filenames.split(" ")]) / (1024**2),
            4,
        )

    def mark_duplicates(
        self, aligned_file: str, out_file: str, metrics_file: str, remove_duplicates: str = "True"
    ) -> str:
        cmd = self.tools.java
        if self.pm.javamem:  # If a memory restriction exists.
            cmd += " -Xmx" + self.pm.javamem
        cmd += " -jar " + self.tools.picard + " MarkDuplicates"
        cmd += " INPUT=" + aligned_file
        cmd += " OUTPUT=" + out_file
        cmd += " METRICS_FILE=" + metrics_file
        cmd += " REMOVE_DUPLICATES=" + remove_duplicates
        return cmd

    def bam2fastq(
        self,
        input_bam: str,
        output_fastq: str,
        output_fastq2: str | None = None,
        unpaired_fastq: str | None = None,
    ) -> str:
        """Build command to convert BAM to FASTQ via Picard SamToFastq."""
        self._ensure_folders(output_fastq, output_fastq2, unpaired_fastq)
        cmd = self.tools.java + " -Xmx" + self.pm.javamem
        cmd += " -jar " + self.tools.picard + " SamToFastq"
        cmd += " INPUT={0}".format(input_bam)
        cmd += " FASTQ={0}".format(output_fastq)
        if output_fastq2 is not None and unpaired_fastq is not None:
            cmd += " SECOND_END_FASTQ={0}".format(output_fastq2)
            cmd += " UNPAIRED_FASTQ={0}".format(unpaired_fastq)
        return cmd

    def bam_to_fastq(self, bam_file: str, out_fastq_pre: str, paired_end: bool) -> str:
        """Build Picard SamToFastq command for BAM to FASTQ conversion."""
        self.make_sure_path_exists(os.path.dirname(out_fastq_pre))
        cmd = self.tools.java + " -Xmx" + self.pm.javamem
        cmd += " -jar " + self.tools.picard + " SamToFastq"
        cmd += " I=" + bam_file
        cmd += " F=" + out_fastq_pre + "_R1.fastq"
        if paired_end:
            cmd += " F2=" + out_fastq_pre + "_R2.fastq"
        cmd += " INCLUDE_NON_PF_READS=true"
        cmd += " QUIET=true"
        cmd += " VERBOSITY=ERROR"
        cmd += " VALIDATION_STRINGENCY=SILENT"
        return cmd

    def bam_to_fastq_awk(
        self, bam_file: str, out_fastq_pre: str, paired_end: bool, zipmode: bool = False
    ) -> tuple[str, str, str | None]:
        """Build fast awk-based BAM to FASTQ conversion command.

        Faster than Picard/bedtools but assumes paired reads are properly
        ordered with no singletons.
        """
        self.make_sure_path_exists(os.path.dirname(out_fastq_pre))
        fq1 = out_fastq_pre + "_R1.fastq"
        fq2 = out_fastq_pre + "_R2.fastq"

        if zipmode:
            fq1 = fq1 + ".gz"
            fq2 = fq2 + ".gz"
            fq1_target = ' | "' + self.ziptool + " -c  > " + fq1 + '"'
            fq2_target = ' | "' + self.ziptool + " -c  > " + fq2 + '"'
        else:
            fq1_target = ' > "' + fq1 + '"'
            fq2_target = ' > "' + fq2 + '"'

        if paired_end:
            cmd = self.tools.samtools + " view " + bam_file + " | awk '"
            cmd += r'{ if (NR%2==1) print "@"$1"/1\n"$10"\n+\n"$11' + fq1_target + ";"
            cmd += r' else print "@"$1"/2\n"$10"\n+\n"$11' + fq2_target + "; }"
            cmd += "'"  # end the awk command
        else:
            fq2 = None
            cmd = self.tools.samtools + " view " + bam_file + " | awk '"
            cmd += r'{ print "@"$1"\n"$10"\n+\n"$11' + fq1_target + "; }"
            cmd += "'"
        return cmd, fq1, fq2

    def bam_to_fastq_bedtools(
        self, bam_file: str, out_fastq_pre: str, paired_end: bool
    ) -> tuple[str, str, str | None]:
        """Build bedtools bamtofastq command for BAM to FASTQ conversion."""
        self.make_sure_path_exists(os.path.dirname(out_fastq_pre))
        fq1 = out_fastq_pre + "_R1.fastq"
        fq2 = None
        cmd = self.tools.bedtools + " bamtofastq -i " + bam_file + " -fq " + fq1 + ".fastq"
        if paired_end:
            fq2 = out_fastq_pre + "_R2.fastq"
            cmd += " -fq2 " + fq2

        return cmd, fq1, fq2

    def get_input_ext(self, input_file: str) -> str:
        """Detect input file type: ".bam", ".fastq.gz", or ".fastq"."""
        if input_file.endswith(".bam"):
            input_ext = ".bam"
        elif input_file.endswith(".fastq.gz") or input_file.endswith(".fq.gz"):
            input_ext = ".fastq.gz"
        elif input_file.endswith(".fastq") or input_file.endswith(".fq"):
            input_ext = ".fastq"
        else:
            errmsg = (
                "'{}'; this pipeline can only deal with .bam, .fastq, or .fastq.gz files".format(
                    input_file
                )
            )
            raise UnsupportedFiletypeException(errmsg)
        return input_ext

    def merge_or_link(
        self, input_args: list, raw_folder: str, local_base: str = "sample"
    ) -> str | list[str]:
        """Standardize inputs by linking or merging .bam/.fastq/.fastq.gz files.

        Example:
            local = tk.merge_or_link([["s1_R1.fq.gz"], ["s1_R2.fq.gz"]], "raw/", "sample1")

        For single files, creates a symlink. For multiple files of the same
        type, merges them (cat for fastq, samtools merge for bam).

        Args:
            input_args: List of input file paths or list of lists (R1/R2).
            raw_folder: Directory for the merge/link output.
            local_base: Base name for output file (usually sample name).
        """
        self.make_sure_path_exists(raw_folder)

        if not isinstance(input_args, list):
            raise TypeError(
                "input_args must be a list of file paths, got {t}. "
                "Pass a list even for single files: ['/path/to/file.fastq']".format(
                    t=type(input_args).__name__
                )
            )

        if any(isinstance(i, list) for i in input_args):
            # We have a list of lists. Process each individually.
            local_input_files = list()
            n_input_files = len(list(filter(bool, input_args)))
            print("Number of input file sets: " + str(n_input_files))

            for input_i, input_arg in enumerate(input_args):
                # Count how many non-null items there are in the list;
                # we only append _R1 (etc.) if there are multiple input files.
                if n_input_files > 1:
                    local_base_extended = local_base + "_R" + str(input_i + 1)
                else:
                    local_base_extended = local_base
                if input_arg:
                    out = self.merge_or_link(input_arg, raw_folder, local_base_extended)

                    print("Local input file: '{}'".format(out))
                    # Make sure file exists:
                    if not os.path.isfile(out):
                        print("Not a file: '{}'".format(out))

                    local_input_files.append(out)

            return local_input_files

        else:
            # We have a list of individual arguments. Merge them.

            if len(input_args) == 1:
                # Only one argument in this list. A single input file; we just link
                # it, regardless of file type:
                # Pull the value out of the list
                input_arg = input_args[0]
                input_ext = self.get_input_ext(input_arg)

                # Convert to absolute path
                if not os.path.isabs(input_arg):
                    input_arg = os.path.abspath(input_arg)

                # Link it to into the raw folder
                local_input_abs = os.path.join(raw_folder, local_base + input_ext)
                self.pm.run(
                    "ln -sf " + input_arg + " " + local_input_abs,
                    target=local_input_abs,
                    shell=True,
                )
                # return the local (linked) filename absolute path
                return local_input_abs

            else:
                # Otherwise, there are multiple inputs.
                # If more than 1 input file is given, then these are to be merged
                # if they are in bam format.
                if all([self.get_input_ext(x) == ".bam" for x in input_args]):
                    sample_merged = local_base + ".merged.bam"
                    output_merge = os.path.join(raw_folder, sample_merged)
                    cmd = self.merge_bams_samtools(input_args, output_merge)
                    self.pm.debug("cmd: {}".format(cmd))
                    self.pm.run(cmd, output_merge)
                    self.validate_bam(output_merge)
                    self.pm.run(cmd, output_merge, nofail=True)
                    return output_merge

                # if multiple fastq
                if all([self.get_input_ext(x) == ".fastq.gz" for x in input_args]):
                    sample_merged_gz = local_base + ".merged.fastq.gz"
                    output_merge_gz = os.path.join(raw_folder, sample_merged_gz)
                    # cmd1 = self.ziptool + "-d -c " + " ".join(input_args) + " > " + output_merge
                    # cmd2 = self.ziptool + " " + output_merge
                    # self.pm.run([cmd1, cmd2], output_merge_gz)
                    # you can save yourself the decompression/recompression:
                    cmd = "cat " + " ".join(input_args) + " > " + output_merge_gz
                    self.pm.run(cmd, output_merge_gz)
                    return output_merge_gz

                if all([self.get_input_ext(x) == ".fastq" for x in input_args]):
                    sample_merged = local_base + ".merged.fastq"
                    output_merge = os.path.join(raw_folder, sample_merged)
                    cmd = "cat " + " ".join(input_args) + " > " + output_merge
                    self.pm.run(cmd, output_merge)
                    return output_merge

                # At this point, we don't recognize the input file types or they
                # do not match.
                raise NotImplementedError(
                    "Cannot merge input files of different types. All input files must be "
                    "either BAM (.bam) or FASTQ (.fastq/.fq/.fastq.gz/.fq.gz). "
                    "Received mixed types. Check your input file list."
                )

    def input_to_fastq(
        self,
        input_file: str | list[str],
        sample_name: str,
        paired_end: bool,
        fastq_folder: str,
        output_file: str | None = None,
        multiclass: bool = False,
        zipmode: bool = False,
    ) -> list:
        """Build command to convert any input (.bam/.fastq.gz/.fastq) to FASTQ.

        Args:
            input_file: Path(s) to input file(s).
            sample_name: Sample name for output file naming.
            paired_end: Whether data is paired-end.
            fastq_folder: Directory for output FASTQ files.
            output_file: Explicit output path (auto-derived if None).
            multiclass: Internal flag for recursive R1/R2 handling.
            zipmode: Output as .fastq.gz instead of .fastq.

        Returns:
            List of [command, fastq_prefix, output_file].
        """

        fastq_prefix = os.path.join(fastq_folder, sample_name)
        self.make_sure_path_exists(fastq_folder)

        # this expects a list; if it gets a string, wrap it in a list.
        if not isinstance(input_file, list):
            input_file = [input_file]

        # If multiple files were provided, recurse on each file individually
        if len(input_file) > 1:
            cmd = []
            output_file = []
            for in_i, in_arg in enumerate(input_file):
                output = fastq_prefix + "_R" + str(in_i + 1) + ".fastq"
                result_cmd, uf, result_file = self.input_to_fastq(
                    in_arg,
                    sample_name,
                    paired_end,
                    fastq_folder,
                    output,
                    multiclass=True,
                    zipmode=zipmode,
                )
                cmd.append(result_cmd)
                output_file.append(result_file)

        else:
            # There was only 1 input class.
            # Convert back into a string
            input_file = input_file[0]
            if not output_file:
                output_file = fastq_prefix + "_R1.fastq"
            if zipmode:
                output_file = output_file + ".gz"

            input_ext = self.get_input_ext(input_file)  # handles .fq or .fastq

            if input_ext == ".bam":
                print("Found .bam file")
                # cmd = self.bam_to_fastq(input_file, fastq_prefix, paired_end)
                cmd, fq1, fq2 = self.bam_to_fastq_awk(
                    input_file, fastq_prefix, paired_end, zipmode
                )
                # pm.run(cmd, output_file, follow=check_fastq)
                if fq2:
                    output_file = [fq1, fq2]
                else:
                    output_file = fq1
            elif input_ext == ".fastq.gz":
                print("Found .fastq.gz file")
                if paired_end and not multiclass:
                    if zipmode:
                        raise NotImplementedError("Can't use zipmode on interleaved fastq data.")
                    # For paired-end reads in one fastq file, we must split the
                    # file into 2. The pipeline author will need to include this
                    # python script in the scripts directory.
                    # TODO: make this self-contained in pypiper. This is a rare
                    # use case these days, as fastq files are almost never
                    # interleaved anymore.
                    script_path = os.path.join(self.tools.scripts_dir, "fastq_split.py")
                    cmd = self.tools.python + " -u " + script_path
                    cmd += " -i " + input_file
                    cmd += " -o " + fastq_prefix
                    # Must also return the set of output files
                    output_file = [
                        fastq_prefix + "_R1.fastq",
                        fastq_prefix + "_R2.fastq",
                    ]
                else:
                    if zipmode:
                        # we do nothing!
                        cmd = "ln -sf " + input_file + " " + output_file
                        print("Found .fq.gz file; no conversion necessary")
                    else:
                        # For single-end reads, we just unzip the fastq.gz file.
                        # or, paired-end reads that were already split.
                        cmd = self.ziptool + " -d -c " + input_file + " > " + output_file
                        # a non-shell version
                        # cmd1 = "gunzip --force " + input_file
                        # cmd2 = "mv " + os.path.splitext(input_file)[0] + " " + output_file
                        # cmd = [cmd1, cmd2]
            elif input_ext == ".fastq":
                if zipmode:
                    cmd = self.ziptool + " -c " + input_file + " > " + output_file
                else:
                    cmd = "ln -sf " + input_file + " " + output_file
                    print("Found .fastq file; no conversion necessary")

        return [cmd, fastq_prefix, output_file]

    def check_fastq(
        self, input_files: str | list[str], output_files: str | list[str], paired_end: bool
    ) -> Callable:
        """Return a follow function that validates FASTQ conversion read counts.

        Example:
            cmd, prefix, out = tk.input_to_fastq(bam, name, pe, folder)
            pm.run(cmd, out, follow=tk.check_fastq(bam, out, pe))

        Args:
            input_files: Original input file(s) before conversion.
            output_files: FASTQ output file(s) from conversion.
            paired_end: Whether data is paired-end.

        Returns:
            Callable that compares read counts and reports stats.
        """

        # Define a temporary function which we will return, to be called by the
        # pipeline.
        # Must define default parameters here based on the parameters passed in. This locks
        # these values in place, so that the variables will be defined when this function
        # is called without parameters as a follow function by pm.run.

        # This is AFTER merge, so if there are multiple files it means the
        # files were split into read1/read2; therefore I must divide by number
        # of files for final reads.
        def temp_func(input_files=input_files, output_files=output_files, paired_end=paired_end):
            if not isinstance(input_files, list):
                input_files = [input_files]
            if not isinstance(output_files, list):
                output_files = [output_files]

            n_input_files = len(list(filter(bool, input_files)))
            n_output_files = len(list(filter(bool, output_files)))

            total_reads = sum(
                [int(self.count_reads(input_file, paired_end)) for input_file in input_files]
            )
            raw_reads = int(total_reads / n_input_files)
            self.pm.report_result("Raw_reads", str(raw_reads))

            total_fastq_reads = sum(
                [int(self.count_reads(output_file, paired_end)) for output_file in output_files]
            )
            fastq_reads = int(total_fastq_reads / n_output_files)

            self.pm.report_result("Fastq_reads", fastq_reads)
            input_ext = self.get_input_ext(input_files[0])
            # We can only assess pass filter reads in bam files with flags.
            if input_ext == ".bam":
                num_failed_filter = sum(
                    [int(self.count_fail_reads(f, paired_end)) for f in input_files]
                )
                pf_reads = int(raw_reads) - num_failed_filter
                self.pm.report_result("PF_reads", str(pf_reads))
            if fastq_reads != int(raw_reads):
                raise Exception(
                    "Fastq conversion error? Number of input reads doesn't number of output reads."
                )

            return fastq_reads

        return temp_func

    def check_trim(
        self,
        trimmed_fastq: str,
        paired_end: bool,
        trimmed_fastq_R2: str | None = None,
        fastqc_folder: str | None = None,
    ) -> Callable:
        """Return a follow function that counts trimmed reads and optionally runs FastQC.

        Args:
            trimmed_fastq: Path to trimmed reads file.
            paired_end: Whether data is paired-end.
            trimmed_fastq_R2: Path to R2 trimmed file for paired-end.
            fastqc_folder: If set, run FastQC and place output here.

        Returns:
            Callable for use as pm.run() follow function.
        """

        def temp_func():
            print("Evaluating read trimming")

            if paired_end and not trimmed_fastq_R2:
                print("WARNING: specified paired-end but no R2 file")

            n_trim = float(self.count_reads(trimmed_fastq, paired_end))
            self.pm.report_result("Trimmed_reads", int(n_trim))
            try:
                rr = float(self.pm.get_stat("Raw_reads"))
            except Exception:
                print("Can't calculate trim loss rate without raw read result.")
            else:
                self.pm.report_result("Trim_loss_rate", round((rr - n_trim) * 100 / rr, 2))

            # Also run a fastqc (if installed/requested)
            if fastqc_folder:
                if fastqc_folder and os.path.isabs(fastqc_folder):
                    self.make_sure_path_exists(fastqc_folder)
                cmd = self.fastqc(trimmed_fastq, fastqc_folder)
                self.pm.run(cmd, lock_name="trimmed_fastqc", nofail=True)
                fname, ext = os.path.splitext(os.path.basename(trimmed_fastq))
                fastqc_html = os.path.join(fastqc_folder, fname + "_fastqc.html")
                self.pm.report_result(
                    "FastQC_report_R1", {"path": fastqc_html, "title": "FastQC report R1"}
                )

                if paired_end and trimmed_fastq_R2:
                    cmd = self.fastqc(trimmed_fastq_R2, fastqc_folder)
                    self.pm.run(cmd, lock_name="trimmed_fastqc_R2", nofail=True)
                    fname, ext = os.path.splitext(os.path.basename(trimmed_fastq_R2))
                    fastqc_html = os.path.join(fastqc_folder, fname + "_fastqc.html")
                    self.pm.report_result(
                        "FastQC_report_R2", {"path": fastqc_html, "title": "FastQC report R2"}
                    )

        return temp_func

    def validate_bam(self, input_bam: str) -> str:
        """Build Picard ValidateSamFile command."""
        cmd = self.tools.java + " -Xmx" + self.pm.javamem
        cmd += " -jar " + self.tools.picard + " ValidateSamFile"
        cmd += " INPUT=" + input_bam
        return cmd

    def merge_bams(
        self,
        input_bams: list[str],
        merged_bam: str,
        in_sorted: bool | str = "TRUE",
        tmp_dir: str | None = None,
    ) -> str | int:
        """Build Picard MergeSamFiles command to combine BAM files."""
        if not len(input_bams) > 1:
            print("No merge required")
            return 0

        outdir, _ = os.path.split(merged_bam)
        if outdir and not os.path.exists(outdir):
            print("Creating path to merge file's folder: '{}'".format(outdir))
            os.makedirs(outdir)

        # Handle more intuitive boolean argument.
        if in_sorted in [False, True]:
            in_sorted = "TRUE" if in_sorted else "FALSE"

        input_string = " INPUT=" + " INPUT=".join(input_bams)
        cmd = self.tools.java + " -Xmx" + self.pm.javamem
        cmd += " -jar " + self.tools.picard + " MergeSamFiles"
        cmd += input_string
        cmd += " OUTPUT=" + merged_bam
        sort_order = "coordinate" if str(in_sorted).upper() == "TRUE" else "unsorted"
        cmd += " ASSUME_SORT_ORDER=" + sort_order
        cmd += " CREATE_INDEX=TRUE"
        cmd += " VALIDATION_STRINGENCY=SILENT"
        if tmp_dir:
            cmd += " TMP_DIR=" + tmp_dir

        return cmd

    def merge_bams_samtools(self, input_bams: list[str], merged_bam: str) -> str:
        cmd = self.tools.samtools + " merge -f "
        cmd += " -@ " + str(self.pm.cores)
        cmd += " " + merged_bam + " "
        cmd += " ".join(input_bams)
        return cmd

    def merge_fastq(
        self, inputs: list[str], output: str, run: bool = False, remove_inputs: bool = False
    ) -> str | None:
        """Merge multiple FASTQ files into one via cat."""
        if remove_inputs and not run:
            raise ValueError("Can't delete files if command isn't run")
        cmd = "cat {} > {}".format(" ".join(inputs), output)
        if run:
            subprocess.check_call(cmd.split(), shell=True)
            if remove_inputs:
                cmd = "rm {}".format(" ".join(inputs))
                subprocess.check_call(cmd.split(), shell=True)
        else:
            return cmd

    def count_lines(self, file_name: str) -> str:
        """Count lines in a file using wc -l."""
        x = subprocess.check_output(
            "wc -l " + file_name + " | sed -E 's/^[[:space:]]+//' | cut -f1 -d' '",
            shell=True,
        )
        return x.decode().strip()

    def count_lines_zip(self, file_name: str) -> str:
        """Count lines in a gzipped file using zcat | wc -l."""
        x = subprocess.check_output(
            self.ziptool
            + " -d -c "
            + file_name
            + " | wc -l | sed -E 's/^[[:space:]]+//' | cut -f1 -d' '",
            shell=True,
        )
        return x.decode().strip()

    def get_chrs_from_bam(self, file_name: str) -> list[str]:
        """Extract chromosome names from a BAM file header via samtools."""
        x = subprocess.check_output(
            self.tools.samtools
            + " view -H "
            + file_name
            + " | grep '^@SQ' | cut -f2| sed s'/SN://'",
            shell=True,
        )
        # Chromosomes will be separated by newlines; split into list to return
        return x.decode().split()

    ###################################
    # Read counting functions
    ###################################
    # In these functions, A paired-end read, with 2 sequences, counts as a two reads

    def count_unique_reads(self, file_name: str, paired_end: bool) -> int:
        """Count unique reads (by name) in a BAM/SAM file. Paired-end counts as 2."""
        if file_name.endswith("sam"):
            param = "-S"
        if file_name.endswith("bam"):
            param = ""
        if paired_end:
            r1 = self.samtools_view(
                file_name,
                param=param + " -f64",
                postpend=" | cut -f1 | sort -k1,1 -u | wc -l | sed -E 's/^[[:space:]]+//'",
            )
            r2 = self.samtools_view(
                file_name,
                param=param + " -f128",
                postpend=" | cut -f1 | sort -k1,1 -u | wc -l | sed -E 's/^[[:space:]]+//'",
            )
        else:
            r1 = self.samtools_view(
                file_name,
                param=param + "",
                postpend=" | cut -f1 | sort -k1,1 -u | wc -l | sed -E 's/^[[:space:]]+//'",
            )
            r2 = 0
        return int(r1) + int(r2)

    def count_unique_mapped_reads(self, file_name: str, paired_end: bool) -> int:
        """Count mapped reads (by name, deduplicated) in a BAM/SAM file."""

        _, ext = os.path.splitext(file_name)
        ext = ext.lower()

        if ext == ".sam":
            param = "-S -F4"
        elif ext == ".bam":
            param = "-F4"
        else:
            raise ValueError(
                "Expected a SAM or BAM file (extension .sam or .bam), "
                "got: '{file}'. Check the file path and extension.".format(file=file_name)
            )

        if paired_end:
            r1 = self.samtools_view(
                file_name,
                param=param + " -f64",
                postpend=" | cut -f1 | sort -k1,1 -u | wc -l | sed -E 's/^[[:space:]]+//'",
            )
            r2 = self.samtools_view(
                file_name,
                param=param + " -f128",
                postpend=" | cut -f1 | sort -k1,1 -u | wc -l | sed -E 's/^[[:space:]]+//'",
            )
        else:
            r1 = self.samtools_view(
                file_name,
                param=param + "",
                postpend=" | cut -f1 | sort -k1,1 -u | wc -l | sed -E 's/^[[:space:]]+//'",
            )
            r2 = 0

        return int(r1) + int(r2)

    def count_flag_reads(self, file_name: str, flag: int | str, paired_end: bool) -> str:
        """Count reads with a specific SAM flag value."""

        param = " -c -f" + str(flag)
        if file_name.endswith("sam"):
            param += " -S"
        return self.samtools_view(file_name, param=param)

    def count_multimapping_reads(self, file_name: str, paired_end: bool) -> int:
        """Count reads flagged as multimapping (SAM flag 256)."""
        return int(self.count_flag_reads(file_name, 256, paired_end))

    def count_uniquelymapping_reads(self, file_name: str, paired_end: bool) -> str:
        """Count reads that mapped to a unique position (exclude flag 256)."""
        param = " -c -F256"
        if file_name.endswith("sam"):
            param += " -S"
        return self.samtools_view(file_name, param=param)

    def count_fail_reads(self, file_name: str, paired_end: bool) -> int:
        """Count reads that failed platform/vendor quality checks (SAM flag 512)."""
        return int(self.count_flag_reads(file_name, 512, paired_end))

    def samtools_view(self, file_name: str, param: str, postpend: str = "") -> str:
        """Run samtools view with given parameters and optional post-processing pipe."""
        cmd = "{} view {} {} {}".format(self.tools.samtools, param, file_name, postpend)
        # in python 3, check_output returns a byte string which causes issues.
        # with python 3.6 we could use argument: "encoding='UTF-8'""
        return subprocess.check_output(cmd, shell=True).decode().strip()

    def count_reads(self, file_name: str, paired_end: bool) -> str | int | float:
        """Count reads in a BAM/SAM/FASTQ file.

        Paired-end reads count as 2. Assumes paired-end FASTQs are split
        into separate R1/R2 files (divides line count by 2 instead of 4).

        Args:
            file_name: Path to BAM/SAM/FASTQ file.
            paired_end: Whether the file contains paired-end reads.
        """

        _, ext = os.path.splitext(file_name)
        if not (is_sam_or_bam(file_name) or is_fastq(file_name)):
            # TODO: make this an exception and force caller to handle that
            # rather than relying on knowledge of possibility of negative value.
            return -1

        if is_sam_or_bam(file_name):
            param_text = "-c" if ext == ".bam" else "-c -S"
            return self.samtools_view(file_name, param=param_text)
        else:
            num_lines = (
                self.count_lines_zip(file_name)
                if is_gzipped_fastq(file_name)
                else self.count_lines(file_name)
            )
            divisor = 2 if paired_end else 4
            return int(num_lines) / divisor

    def count_concordant(self, aligned_bam: str) -> str:
        """Count reads aligned concordantly exactly once (YT:Z:CP flag)."""
        cmd = self.tools.samtools + " view " + aligned_bam + " | "
        cmd += "grep 'YT:Z:CP'" + " | uniq -u | wc -l | sed -E 's/^[[:space:]]+//'"

        return subprocess.check_output(cmd, shell=True).decode().strip()

    def count_mapped_reads(self, file_name: str, paired_end: bool) -> str | int:
        """Count mapped reads in a BAM/SAM file (excludes unmapped, flag -F4)."""
        if file_name.endswith("bam"):
            return self.samtools_view(file_name, param="-c -F4")
        if file_name.endswith("sam"):
            return self.samtools_view(file_name, param="-c -F4 -S")
        return -1

    def sam_conversions(self, sam_file: str, depth: bool = True) -> str:
        """Build command to convert SAM to sorted/indexed BAM (optionally with depth).

        Args:
            sam_file: Path to SAM file.
            depth: Also calculate per-position coverage.
        """
        cmd = (
            self.tools.samtools
            + " view -bS "
            + sam_file
            + " > "
            + sam_file.replace(".sam", ".bam")
            + "\n"
        )
        cmd += (
            self.tools.samtools
            + " sort "
            + sam_file.replace(".sam", ".bam")
            + " -o "
            + sam_file.replace(".sam", "_sorted.bam")
            + "\n"
        )
        cmd += self.tools.samtools + " index " + sam_file.replace(".sam", "_sorted.bam") + "\n"
        if depth:
            cmd += (
                self.tools.samtools
                + " depth "
                + sam_file.replace(".sam", "_sorted.bam")
                + " > "
                + sam_file.replace(".sam", "_sorted.depth")
                + "\n"
            )
        return cmd

    def bam_conversions(self, bam_file: str, depth: bool = True) -> str:
        """Build command to sort and index a BAM file (optionally with depth).

        Args:
            bam_file: Path to BAM file.
            depth: Also calculate per-position coverage.
        """
        cmd = (
            self.tools.samtools
            + " view -h "
            + bam_file
            + " > "
            + bam_file.replace(".bam", ".sam")
            + "\n"
        )
        cmd += (
            self.tools.samtools
            + " sort "
            + bam_file
            + " -o "
            + bam_file.replace(".bam", "_sorted.bam")
            + "\n"
        )
        cmd += self.tools.samtools + " index " + bam_file.replace(".bam", "_sorted.bam") + "\n"
        if depth:
            cmd += (
                self.tools.samtools
                + " depth "
                + bam_file.replace(".bam", "_sorted.bam")
                + " > "
                + bam_file.replace(".bam", "_sorted.depth")
                + "\n"
            )
        return cmd

    def fastqc(self, file: str, output_dir: str) -> str:
        """Build FastQC command for a reads file."""
        # You can find the fastqc help with fastqc --help
        try:
            pm = self.pm
        except AttributeError:
            # Do nothing, this is just for path construction.
            pass
        else:
            if not os.path.isabs(output_dir) and pm is not None:
                output_dir = os.path.join(pm.outfolder, output_dir)
        self.make_sure_path_exists(output_dir)
        return "{} --noextract --outdir {} {}".format(self.tools.fastqc, output_dir, file)

    def fastqc_rename(self, input_bam: str, output_dir: str, sample_name: str) -> list[str]:
        """Build commands to run FastQC and rename output by sample name."""
        cmds = list()
        initial = os.path.splitext(os.path.basename(input_bam))[0]
        cmd1 = self.fastqc(input_bam, output_dir)
        cmds.append(cmd1)
        cmd2 = "if [[ ! -s {1}_fastqc.html ]]; then mv {0}_fastqc.html {1}_fastqc.html; mv {0}_fastqc.zip {1}_fastqc.zip; fi".format(
            os.path.join(output_dir, initial), os.path.join(output_dir, sample_name)
        )
        cmds.append(cmd2)
        return cmds

    def samtools_index(self, bam_file: str) -> str:
        """Index a bam file."""
        cmd = self.tools.samtools + " index {0}".format(bam_file)
        return cmd

    def slurm_header(
        self,
        job_name: str,
        output: str,
        queue: str = "shortq",
        n_tasks: int = 1,
        time: str = "10:00:00",
        cpus_per_task: int = 8,
        mem_per_cpu: int = 2000,
        nodes: int = 1,
        user_mail: str = "",
        mail_type: str = "end",
    ) -> str:
        cmd = """       #!/bin/bash
        #SBATCH --partition={0}
        #SBATCH --ntasks={1}
        #SBATCH --time={2}

        #SBATCH --cpus-per-task={3}
        #SBATCH --mem-per-cpu={4}
        #SBATCH --nodes={5}

        #SBATCH --job-name={6}
        #SBATCH --output={7}

        #SBATCH --mail-type={8}
        #SBATCH --mail-user={9}

        # Start running the job
        hostname
        date

        """.format(
            queue,
            n_tasks,
            time,
            cpus_per_task,
            mem_per_cpu,
            nodes,
            job_name,
            output,
            mail_type,
            user_mail,
        )

        return cmd

    def slurm_footer(self) -> str:
        return "     date"

    def slurm_submit_job(self, job_file: str) -> int:
        return os.system("sbatch %s" % job_file)

    def remove_file(self, file_name: str) -> str:
        return "rm {0}".format(file_name)

    def move_file(self, old: str, new: str) -> str:
        return "mv {0} {1}".format(old, new)

    def preseq_curve(self, bam_file: str, output_prefix: str) -> str:
        return """
        preseq c_curve -B -P -o {0}.yield.txt {1}
        """.format(output_prefix, bam_file)

    def preseq_extrapolate(self, bam_file: str, output_prefix: str) -> str:
        return """
        preseq lc_extrap -v -B -P -e 1e+9 -o {0}.future_yield.txt {1}
        """.format(output_prefix, bam_file)

    def preseq_coverage(self, bam_file: str, output_prefix: str) -> str:
        return """
        preseq gc_extrap -o {0}.future_coverage.txt {1}
        """.format(output_prefix, bam_file)

    def trimmomatic(
        self,
        input_fastq1: str,
        output_fastq1: str,
        cpus: int | str,
        adapters: str,
        log: str,
        input_fastq2: str | None = None,
        output_fastq1_unpaired: str | None = None,
        output_fastq2: str | None = None,
        output_fastq2_unpaired: str | None = None,
    ) -> str:
        PE = False if input_fastq2 is None else True
        pe = "PE" if PE else "SE"
        cmd = self.tools.java + " -Xmx" + self.pm.javamem
        cmd += " -jar " + self.tools.trimmomatic
        cmd += " {0} -threads {1} -trimlog {2} {3}".format(pe, cpus, log, input_fastq1)
        if PE:
            cmd += " {0}".format(input_fastq2)
        cmd += " {0}".format(output_fastq1)
        if PE:
            cmd += " {0} {1} {2}".format(
                output_fastq1_unpaired, output_fastq2, output_fastq2_unpaired
            )
        cmd += " ILLUMINACLIP:{0}:1:40:15:8:true".format(adapters)
        cmd += " LEADING:3 TRAILING:3"
        cmd += " SLIDINGWINDOW:4:10"
        cmd += " MINLEN:36"
        return cmd

    def skewer(
        self,
        input_fastq1: str,
        output_prefix: str,
        output_fastq1: str,
        log: str,
        cpus: int | str,
        adapters: str,
        input_fastq2: str | None = None,
        output_fastq2: str | None = None,
    ) -> list[str]:
        """Build skewer adapter-trimming commands with file renaming."""

        pe = input_fastq2 is not None
        mode = "pe" if pe else "any"
        cmds = list()
        cmd1 = self.tools.skewer + " --quiet"
        cmd1 += " -f sanger"
        cmd1 += " -t {0}".format(cpus)
        cmd1 += " -m {0}".format(mode)
        cmd1 += " -x {0}".format(adapters)
        cmd1 += " -o {0}".format(output_prefix)
        cmd1 += " {0}".format(input_fastq1)
        if input_fastq2 is None:
            cmds.append(cmd1)
        else:
            cmd1 += " {0}".format(input_fastq2)
            cmds.append(cmd1)
        if input_fastq2 is None:
            cmd2 = "mv {0} {1}".format(output_prefix + "-trimmed.fastq", output_fastq1)
            cmds.append(cmd2)
        else:
            cmd2 = "mv {0} {1}".format(output_prefix + "-trimmed-pair1.fastq", output_fastq1)
            cmds.append(cmd2)
            cmd3 = "mv {0} {1}".format(output_prefix + "-trimmed-pair2.fastq", output_fastq2)
            cmds.append(cmd3)
        cmd4 = "mv {0} {1}".format(output_prefix + "-trimmed.log", log)
        cmds.append(cmd4)
        return cmds

    def bowtie2_map(
        self,
        input_fastq1: str,
        output_bam: str,
        log: str,
        metrics: str,
        genome_index: str,
        max_insert: int | str,
        cpus: int | str,
        input_fastq2: str | None = None,
    ) -> str:
        # Admits 2000bp-long fragments (--maxins option)
        cmd = self.tools.bowtie2 + " --very-sensitive --no-discordant -p {0}".format(cpus)
        cmd += " -x {0}".format(genome_index)
        cmd += " --met-file {0}".format(metrics)
        if input_fastq2 is None:
            cmd += " {0} ".format(input_fastq1)
        else:
            cmd += " --maxins {0}".format(max_insert)
            cmd += " -1 {0}".format(input_fastq1)
            cmd += " -2 {0}".format(input_fastq2)
        cmd += " 2> {0} | samtools view -S -b - | samtools sort -o {1} -".format(log, output_bam)
        return cmd

    def topHat_map(
        self, input_fastq: str, output_dir: str, genome: str, transcriptome: str, cpus: int | str
    ) -> str:
        # TODO:
        # Allow paired input
        cmd = (
            self.tools.tophat
            + " --GTF {0} --b2-L 15 --library-type fr-unstranded --mate-inner-dist 120".format(
                transcriptome
            )
        )
        cmd += " --max-multihits 100 --no-coverage-search"
        cmd += " --num-threads {0} --output-dir {1} {2} {3}".format(
            cpus, output_dir, genome, input_fastq
        )
        return cmd

    def picard_mark_duplicates(
        self, input_bam: str, output_bam: str, metrics_file: str, temp_dir: str = "."
    ) -> list[str]:
        transient_file = re.sub(r"\.bam$", "", output_bam) + ".dups.nosort.bam"
        output_bam = re.sub(r"\.bam$", "", output_bam)
        cmd1 = self.tools.java + " -Xmx" + self.pm.javamem
        cmd1 += " -jar  `which MarkDuplicates.jar`"
        cmd1 += " INPUT={0}".format(input_bam)
        cmd1 += " OUTPUT={0}".format(transient_file)
        cmd1 += " METRICS_FILE={0}".format(metrics_file)
        cmd1 += " VALIDATION_STRINGENCY=LENIENT"
        cmd1 += " TMP_DIR={0}".format(temp_dir)
        # Sort bam file with marked duplicates
        cmd2 = self.tools.samtools + " sort {0} {1}".format(transient_file, output_bam)
        # Remove transient file
        cmd3 = "if [[ -s {0} ]]; then rm {0}; fi".format(transient_file)
        return [cmd1, cmd2, cmd3]

    def sambamba_remove_duplicates(self, input_bam: str, output_bam: str, cpus: int = 16) -> str:
        cmd = self.tools.sambamba + " markdup -t {0} -r {1} {2}".format(
            cpus, input_bam, output_bam
        )
        return cmd

    def get_mitochondrial_reads(self, bam_file: str, output: str, cpus: int = 4) -> list[str]:
        """ """
        tmp_bam = bam_file + "tmp_rmMe"
        cmd1 = self.tools.sambamba + " index -t {0} {1}".format(cpus, bam_file)
        cmd2 = (
            self.tools.sambamba
            + " slice {0} chrM | {1} markdup -t 4 /dev/stdin {2} 2> {3}".format(
                bam_file, self.tools.sambamba, tmp_bam, output
            )
        )
        cmd3 = "rm {}".format(tmp_bam)
        return [cmd1, cmd2, cmd3]

    def filter_reads(
        self,
        input_bam: str,
        output_bam: str,
        metrics_file: str,
        paired: bool = False,
        cpus: int = 16,
        Q: int = 30,
    ) -> list[str]:
        """Build commands to dedup, quality-filter, and remove multimappers."""
        nodups = re.sub(r"\.bam$", "", output_bam) + ".nodups.nofilter.bam"
        cmd1 = (
            self.tools.sambamba
            + " markdup -t {0} -r --compression-level=0 {1} {2} 2> {3}".format(
                cpus, input_bam, nodups, metrics_file
            )
        )
        cmd2 = self.tools.sambamba + " view -t {0} -f bam --valid".format(cpus)
        if paired:
            cmd2 += ' -F "not (unmapped or mate_is_unmapped) and proper_pair'
        else:
            cmd2 += ' -F "not unmapped'
        cmd2 += (
            ' and not (secondary_alignment or supplementary) and mapping_quality >= {0}"'.format(Q)
        )
        cmd2 += " {0} |".format(nodups)
        cmd2 += self.tools.sambamba + " sort -t {0} /dev/stdin -o {1}".format(cpus, output_bam)
        cmd3 = "if [[ -s {0} ]]; then rm {0}; fi".format(nodups)
        cmd4 = "if [[ -s {0} ]]; then rm {0}; fi".format(nodups + ".bai")
        return [cmd1, cmd2, cmd3, cmd4]

    def shift_reads(self, input_bam: str, genome: str, output_bam: str) -> str:
        # output_bam = re.sub("\.bam$", "", output_bam)
        cmd = self.tools.samtools + " view -h {0} |".format(input_bam)
        cmd += " shift_reads.py {0} |".format(genome)
        cmd += " " + self.tools.samtools + " view -S -b - |"
        cmd += " " + self.tools.samtools + " sort -o {0} -".format(output_bam)
        return cmd

    def sort_index_bam(self, input_bam: str, output_bam: str) -> list[str]:
        tmp_bam = re.sub(r"\.bam", ".sorted", input_bam)
        cmd1 = self.tools.samtools + " sort {0} {1}".format(input_bam, tmp_bam)
        cmd2 = "mv {0}.bam {1}".format(tmp_bam, output_bam)
        cmd3 = self.tools.samtools + " index {0}".format(output_bam)
        return [cmd1, cmd2, cmd3]

    def index_bam(self, input_bam: str) -> str:
        cmd = self.tools.samtools + " index {0}".format(input_bam)
        return cmd

    def run_spp(self, input_bam: str, output: str, plot: str, cpus: int) -> str:
        """
        Run the SPP read peak analysis tool.

        Args:
            input_bam (str): Path to reads file
            output (str): Path to output file
            plot (str): Path to plot file
            cpus (int): Number of processors to use

        Returns:
            str: Command with which to run SPP
        """
        base = "{} {} -rf -savp".format(self.tools.Rscript, self.tools.spp)
        cmd = base + " -savp={} -s=0:5:500 -c={} -out={} -p={}".format(
            plot, input_bam, output, cpus
        )
        return cmd

    def get_fragment_sizes(self, bam_file: str) -> Any:
        try:
            import numpy as np
            import pysam
        except Exception:
            return
        frag_sizes = list()
        bam = pysam.Samfile(bam_file, "rb")
        for read in bam:
            if bam.getrname(read.tid) != "chrM" and read.tlen < 1500:
                frag_sizes.append(read.tlen)
        bam.close()
        return np.array(frag_sizes)

    def plot_atacseq_insert_sizes(
        self,
        bam: str,
        plot: str,
        output_csv: str,
        max_insert: int = 1500,
        smallest_insert: int = 30,
    ) -> None:
        """
        Heavy inspiration from here:
        https://github.com/dbrg77/ATAC/blob/master/ATAC_seq_read_length_curve_fitting.ipynb
        """
        try:
            import matplotlib
            import matplotlib.mlab as mlab
            import numpy as np
            import pysam
            from scipy.integrate import simps
            from scipy.optimize import curve_fit

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except Exception:
            print("Necessary Python modules couldn't be loaded.")
            return

        try:
            import seaborn as sns

            sns.set_style("whitegrid")
        except Exception:
            pass

        def get_fragment_sizes(bam, max_insert=1500):
            frag_sizes = list()

            bam = pysam.Samfile(bam, "rb")

            for i, read in enumerate(bam):
                if read.tlen < max_insert:
                    frag_sizes.append(read.tlen)
            bam.close()

            return np.array(frag_sizes)

        def mixture_function(x, *p):
            """
            Mixture function to model four gaussian (nucleosomal)
            and one exponential (nucleosome-free) distributions.
            """
            m1, s1, w1, m2, s2, w2, m3, s3, w3, m4, s4, w4, q, r = p
            nfr = expo(x, 2.9e-02, 2.8e-02)
            nfr[:smallest_insert] = 0

            return (
                mlab.normpdf(x, m1, s1) * w1
                + mlab.normpdf(x, m2, s2) * w2
                + mlab.normpdf(x, m3, s3) * w3
                + mlab.normpdf(x, m4, s4) * w4
                + nfr
            )

        def expo(x, q, r):
            """
            Exponential function.
            """
            return q * np.exp(-r * x)

        # get fragment sizes
        frag_sizes = get_fragment_sizes(bam)

        # bin
        numBins = np.linspace(0, max_insert, max_insert + 1)
        y, scatter_x = np.histogram(frag_sizes, numBins, density=1)
        # get the mid-point of each bin
        x = (scatter_x[:-1] + scatter_x[1:]) / 2

        # Parameters are empirical, need to check
        paramGuess = [
            200,
            50,
            0.7,  # gaussians
            400,
            50,
            0.15,
            600,
            50,
            0.1,
            800,
            55,
            0.045,
            2.9e-02,
            2.8e-02,  # exponential
        ]

        try:
            popt3, pcov3 = curve_fit(
                mixture_function,
                x[smallest_insert:],
                y[smallest_insert:],
                p0=paramGuess,
                maxfev=100000,
            )
        except Exception:
            print("Nucleosomal fit could not be found.")
            return

        m1, s1, w1, m2, s2, w2, m3, s3, w3, m4, s4, w4, q, r = popt3

        # Plot
        plt.figure(figsize=(12, 12))

        # Plot distribution
        plt.hist(frag_sizes, numBins, histtype="step", ec="k", normed=1, alpha=0.5)

        # Plot nucleosomal fits
        plt.plot(x, mlab.normpdf(x, m1, s1) * w1, "r-", lw=1.5, label="1st nucleosome")
        plt.plot(x, mlab.normpdf(x, m2, s2) * w2, "g-", lw=1.5, label="2nd nucleosome")
        plt.plot(x, mlab.normpdf(x, m3, s3) * w3, "b-", lw=1.5, label="3rd nucleosome")
        plt.plot(x, mlab.normpdf(x, m4, s4) * w4, "c-", lw=1.5, label="4th nucleosome")

        # Plot nucleosome-free fit
        nfr = expo(x, 2.9e-02, 2.8e-02)
        nfr[:smallest_insert] = 0
        plt.plot(x, nfr, "k-", lw=1.5, label="nucleosome-free")

        # Plot sum of fits
        ys = mixture_function(x, *popt3)
        plt.plot(x, ys, "k--", lw=3.5, label="fit sum")

        plt.legend()
        plt.xlabel("Fragment size (bp)")
        plt.ylabel("Density")
        plt.savefig(plot, bbox_inches="tight")

        # Integrate curves and get areas under curve
        areas = [
            ["fraction", "area under curve", "max density"],
            ["Nucleosome-free fragments", simps(nfr), max(nfr)],
            [
                "1st nucleosome",
                simps(mlab.normpdf(x, m1, s1) * w1),
                max(mlab.normpdf(x, m1, s1) * w1),
            ],
            [
                "2nd nucleosome",
                simps(mlab.normpdf(x, m2, s2) * w1),
                max(mlab.normpdf(x, m2, s2) * w2),
            ],
            [
                "3rd nucleosome",
                simps(mlab.normpdf(x, m3, s3) * w1),
                max(mlab.normpdf(x, m3, s3) * w3),
            ],
            [
                "4th nucleosome",
                simps(mlab.normpdf(x, m4, s4) * w1),
                max(mlab.normpdf(x, m4, s4) * w4),
            ],
        ]

        try:
            import csv

            with open(output_csv, "w") as f:
                writer = csv.writer(f)
                writer.writerows(areas)
        except Exception:
            pass

    # TODO: parameterize in terms of normalization factor.
    def bam_to_bigwig(
        self,
        input_bam: str,
        output_bigwig: str,
        genome_sizes: str,
        genome: str,
        tagmented: bool = False,
        normalize: bool = False,
        norm_factor: int = 1000,
    ) -> list[str]:
        """
        Convert a BAM file to a bigWig file.

        Args:
            input_bam (str): path to BAM file to convert
            output_bigwig (str): path to which to write file in bigwig format
            genome_sizes (str): path to file with chromosome size information
            genome (str): name of genomic assembly
            tagmented (bool): flag related to read-generating protocol
            normalize (bool): whether to normalize coverage
            norm_factor (int): number of bases to use for normalization

        Returns:
            list[str]: sequence of commands to execute
        """
        # TODO:
        # addjust fragment length dependent on read size and real fragment size
        # (right now it asssumes 50bp reads with 180bp fragments)
        cmds = list()
        transient_file = os.path.abspath(re.sub(r"\.bigWig", "", output_bigwig))
        cmd1 = self.tools.bedtools + " bamtobed -i {0} |".format(input_bam)
        if not tagmented:
            cmd1 += (
                " "
                + self.tools.bedtools
                + " slop -i stdin -g {0} -s -l 0 -r 130 |".format(genome_sizes)
            )
            cmd1 += " fix_bedfile_genome_boundaries.py {0} |".format(genome)
        cmd1 += (
            " "
            + self.tools.genomeCoverageBed
            + " {0}-bg -g {1} -i stdin > {2}.cov".format(
                "-5 " if tagmented else "", genome_sizes, transient_file
            )
        )
        cmds.append(cmd1)
        if normalize:
            cmds.append(
                """awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * {1}; print}}' {0}.cov {0}.cov | sort -k1,1 -k2,2n > {0}.normalized.cov""".format(
                    transient_file, norm_factor
                )
            )
        cmds.append(
            self.tools.bedGraphToBigWig
            + " {0}{1}.cov {2} {3}".format(
                transient_file,
                ".normalized" if normalize else "",
                genome_sizes,
                output_bigwig,
            )
        )
        # remove tmp files
        cmds.append("if [[ -s {0}.cov ]]; then rm {0}.cov; fi".format(transient_file))
        if normalize:
            cmds.append(
                "if [[ -s {0}.normalized.cov ]]; then rm {0}.normalized.cov; fi".format(
                    transient_file
                )
            )
        cmds.append("chmod 755 {0}".format(output_bigwig))
        return cmds

    def add_track_to_hub(
        self, sample_name: str, track_url: str, track_hub: str, colour: str, five_prime: str = ""
    ) -> list[str]:
        cmd1 = """echo "track type=bigWig name='{0} {1}' description='{0} {1}'""".format(
            sample_name, five_prime
        )
        cmd1 += """ height=32 visibility=full maxHeightPixels=32:32:25 bigDataUrl={0} color={1}" >> {2}""".format(
            track_url, colour, track_hub
        )
        cmd2 = "chmod 755 {0}".format(track_hub)
        return [cmd1, cmd2]

    def link_to_track_hub(self, track_hub_url: str, file_name: str, genome: str) -> None:
        import textwrap

        db = "org" if genome == "hg19" else "db"  # different database call for human
        genome = "human" if genome == "hg19" else genome  # change hg19 to human
        html = """
        <html>
            <head>
                <meta http-equiv="refresh" content="0; url=http://genome.ucsc.edu/cgi-bin/hgTracks?"""
        html += """{db}={genome}&hgt.customText={track_hub_url}" />
            </head>
        </html>
        """.format(track_hub_url=track_hub_url, genome=genome, db=db)
        with open(file_name, "w") as handle:
            handle.write(textwrap.dedent(html))

    def htseq_count(self, input_bam: str, gtf: str, output: str) -> list[str]:
        sam = input_bam.replace("bam", "sam")
        cmd1 = "samtools view {0} > {1}".format(input_bam, sam)
        cmd2 = "htseq-count -f sam -t exon -i transcript_id -m union {0} {1} > {2}".format(
            sam, gtf, output
        )
        cmd3 = "rm {0}".format(sam)
        return [cmd1, cmd2, cmd3]

    def kallisto(
        self,
        input_fastq: str,
        output_dir: str,
        output_bam: str,
        transcriptome_index: str,
        cpus: int,
        input_fastq2: str | None = None,
        size: int = 180,
        b: int = 200,
    ) -> list[str]:
        cmd1 = (
            self.tools.kallisto
            + " quant --bias --pseudobam -b {0} -l {1} -i {2} -o {3} -t {4}".format(
                b, size, transcriptome_index, output_dir, cpus
            )
        )
        if input_fastq2 is None:
            cmd1 += " --single {0}".format(input_fastq)
        else:
            cmd1 += " {0} {1}".format(input_fastq, input_fastq2)
        cmd1 += " | " + self.tools.samtools + " view -Sb - > {0}".format(output_bam)
        cmd2 = self.tools.kallisto + " h5dump -o {0} {0}/abundance.h5".format(output_dir)
        return [cmd1, cmd2]

    def genome_wide_coverage(self, input_bam: str, genome_windows: str, output: str) -> str:
        cmd = self.tools.bedtools + " coverage -counts -abam {0} -b {1} > {2}".format(
            input_bam, genome_windows, output
        )
        return cmd

    def calc_frip(self, input_bam: str, input_bed: str, threads: int = 4) -> str:
        """
        Calculate fraction of reads in peaks.

        A file of with a pool of sequencing reads and a file with peak call
        regions define the operation that will be performed. Thread count
        for samtools can be specified as well.

        Args:
            input_bam (str): sequencing reads file
            input_bed (str): file with called peak regions
            threads (int): number of threads samtools may use

        Returns:
            float: fraction of reads in peaks defined in given peaks file
        """
        cmd = self.simple_frip(input_bam, input_bed, threads)
        return subprocess.check_output(cmd.split(" "), shell=True).decode().strip()

    def simple_frip(self, input_bam: str, input_bed: str, threads: int = 4) -> str:
        cmd = "{} view".format(self.tools.samtools)
        cmd += " -@ {} -c -L {}".format(threads, input_bed)
        cmd += " " + input_bam
        return cmd

    def calculate_frip(self, input_bam: str, input_bed: str, output: str, cpus: int = 4) -> str:
        cmd = self.tools.sambamba + " depth region -t {0}".format(cpus)
        cmd += " -L {0}".format(input_bed)
        cmd += " {0}".format(input_bam)
        cmd += " | awk '{{sum+=$5}} END {{print sum}}' > {0}".format(output)
        return cmd

    def macs2_call_peaks(
        self,
        treatment_bams: str | Iterable[str],
        output_dir: str,
        sample_name: str,
        genome: str,
        control_bams: str | Iterable[str] | None = None,
        broad: bool = False,
        paired: bool = False,
        pvalue: float | None = None,
        qvalue: float | None = None,
        include_significance: bool | None = None,
    ) -> str:
        """
        Use MACS2 to call peaks.

        Args:
            treatment_bams (str | Iterable[str]): Paths to files with data to
                regard as treatment.
            output_dir (str): Path to output folder.
            sample_name (str): Name for the sample involved.
            genome (str): Name of the genome assembly to use.
            control_bams (str | Iterable[str]): Paths to files with data to
                regard as control
            broad (bool): Whether to do broad peak calling.
            paired (bool): Whether reads are paired-end
            pvalue (float | NoneType): Statistical significance measure to
                pass as --pvalue to peak calling with MACS
            qvalue (float | NoneType): Statistical significance measure to
                pass as --qvalue to peak calling with MACS
            include_significance (bool | NoneType): Whether to pass a
                statistical significance argument to peak calling with MACS; if
                omitted, this will be True if the peak calling is broad or if
                either p-value or q-value is specified; default significance
                specification is a p-value of 0.001 if a significance is to be
                specified but no value is provided for p-value or q-value.

        Returns:
            str: Command to run.
        """
        sizes = {
            "hg38": 2.7e9,
            "hg19": 2.7e9,
            "mm10": 1.87e9,
            "dr7": 1.412e9,
            "mm9": 1.87e9,
        }

        # Whether to specify to MACS2 a value for statistical significance
        # can be either directly indicated, but if not, it's determined by
        # whether the mark is associated with broad peaks. By default, we
        # specify a significance value to MACS2 for a mark associated with a
        # broad peak.
        if include_significance is None:
            include_significance = broad

        cmd = self.tools.macs2 + " callpeak -t {0}".format(
            treatment_bams if type(treatment_bams) is str else " ".join(treatment_bams)
        )

        if control_bams is not None:
            cmd += " -c {0}".format(
                control_bams if type(control_bams) is str else " ".join(control_bams)
            )

        if paired:
            cmd += " -f BAMPE "

        # Additional settings based on whether the marks is associated with
        # broad peaks
        if broad:
            cmd += " --broad --nomodel --extsize 73"
        else:
            cmd += " --fix-bimodal --extsize 180 --bw 200"

        if include_significance:
            # Allow significance specification via either p- or q-value,
            # giving preference to q-value if both are provided but falling
            # back on a default p-value if neither is provided but inclusion
            # of statistical significance measure is desired.
            if qvalue is not None:
                cmd += " --qvalue {}".format(qvalue)
            else:
                cmd += " --pvalue {}".format(pvalue or 0.00001)
        cmd += " -g {0} -n {1} --outdir {2}".format(sizes[genome], sample_name, output_dir)

        return cmd

    def macs2_call_peaks_atacseq(
        self, treatment_bam: str, output_dir: str, sample_name: str, genome: str
    ) -> str:
        genome_sizes = {
            "hg38": 2.7e9,
            "hg19": 2.7e9,
            "mm10": 1.87e9,
            "dr7": 1.412e9,
            "mm9": 1.87e9,
        }
        cmd = self.tools.macs2 + " callpeak -t {0}".format(treatment_bam)
        cmd += " --nomodel --extsize 147 -g {0} -n {1} --outdir {2}".format(
            genome_sizes[genome], sample_name, output_dir
        )
        return cmd

    def macs2_plot_model(
        self, r_peak_model_file: str, sample_name: str, output_dir: str
    ) -> list[str]:
        # run macs r script
        cmd1 = "{} {}".format(self.tools.Rscript, r_peak_model_file)
        # move output plot to sample dir
        cmd2 = "mv {0}/{1}_model.pdf {2}/{1}_model.pdf".format(
            os.getcwd(), sample_name, output_dir
        )
        return [cmd1, cmd2]

    def spp_call_peaks(
        self,
        treatment_bam: str,
        control_bam: str,
        treatment_name: str,
        control_name: str,
        output_dir: str,
        broad: str | bool,
        cpus: int,
        qvalue: float | None = None,
    ) -> str:
        """
        Build command for R script to call peaks with SPP.

        Args:
            treatment_bam (str): Path to file with data for treatment sample.
            control_bam (str): Path to file with data for control sample.
            treatment_name (str): Name for the treatment sample.
            control_name (str): Name for the control sample.
            output_dir (str): Path to folder for output.
            broad (str | bool): Whether to specify broad peak calling mode.
            cpus (int): Number of cores the script may use.
            qvalue (float): FDR, as decimal value

        Returns:
            str: Command to run.
        """
        broad = "TRUE" if broad else "FALSE"
        cmd = (
            self.tools.Rscript
            + " `which spp_peak_calling.R` {0} {1} {2} {3} {4} {5} {6}".format(
                treatment_bam,
                control_bam,
                treatment_name,
                control_name,
                broad,
                cpus,
                output_dir,
            )
        )
        if qvalue is not None:
            cmd += " {}".format(qvalue)
        return cmd

    def bam_to_bed(self, input_bam: str, output_bed: str) -> str:
        cmd = self.tools.bedtools + " bamtobed -i {0} > {1}".format(input_bam, output_bed)
        return cmd

    def zinba_call_peaks(
        self, treatment_bed: str, control_bed: str, cpus: int, tagmented: bool = False
    ) -> str:
        fragmentLength = 80 if tagmented else 180
        cmd = self.tools.Rscript + " `which zinba.R` -l {0} -t {1} -c {2}".format(
            fragmentLength, treatment_bed, control_bed
        )
        return cmd

    def filter_peaks_mappability(self, peaks: str, alignability: str, filtered_peaks: str) -> str:
        cmd = self.tools.bedtools + " intersect -wa -u -f 1"
        cmd += " -a {0} -b {1} > {2} ".format(peaks, alignability, filtered_peaks)
        return cmd

    def homer_find_motifs(
        self,
        peak_file: str,
        genome: str,
        output_dir: str,
        size: int = 150,
        length: str = "8,10,12,14,16",
        n_motifs: int = 12,
    ) -> str:
        cmd = "findMotifsGenome.pl {0} {1} {2}".format(peak_file, genome, output_dir)
        cmd += " -mask -size {0} -len {1} -S {2}".format(size, length, n_motifs)
        return cmd

    def homer_annotate_pPeaks(
        self, peak_file: str, genome: str, motif_file: str, output_bed: str
    ) -> str:
        cmd = "annotatePeaks.pl {0} {1} -mask -mscore -m {2} |".format(
            peak_file, genome, motif_file
        )
        cmd += "tail -n +2 | cut -f 1,5,22 > {3}"
        return cmd

    def center_peaks_on_motifs(
        self, peak_file: str, genome: str, window_width: int, motif_file: str, output_bed: str
    ) -> str:
        cmd = "annotatePeaks.pl {0} {1} -size {2} -center {3} |".format(
            peak_file, genome, window_width, motif_file
        )
        cmd += " awk -v OFS='\t' '{print $2, $3, $4, $1, $6, $5}' |"
        cmd += (
            """ awk -v OFS='\t' -F '\t' '{ gsub("0", "+", $6) ; gsub("1", "-", $6) ; print }' |"""
        )
        cmd += " fix_bedfile_genome_boundaries.py {0} | sortBed > {1}".format(genome, output_bed)
        return cmd

    def get_read_type(self, bam_file: str, n: int = 10) -> tuple[str, int]:
        """
        Gets the read type (single, paired) and length of bam file.

        Args:
            bam_file (str): Bam file to determine read attributes.
            n (int): Number of lines to read from bam file.

        Returns:
            str, int: tuple of read type and read length
        """

        from collections.abc import Counter

        try:
            p = subprocess.Popen([self.tools.samtools, "view", bam_file], stdout=subprocess.PIPE)
            # Count paired alignments
            paired = 0
            read_length = Counter()
            while n > 0:
                line = p.stdout.next().split("\t")
                flag = int(line[1])
                read_length[len(line[9])] += 1
                if 1 & flag:  # check decimal flag contains 1 (paired)
                    paired += 1
                n -= 1
            p.kill()
        except IOError("Cannot read provided bam file.") as e:
            raise e
        # Get most abundant read read_length
        read_length = sorted(read_length)[-1]
        # If at least half is paired, return True
        if paired > (n / 2.0):
            return "PE", read_length
        else:
            return "SE", read_length

    def parse_bowtie_stats(self, stats_file: str) -> dict:
        """
        Parses Bowtie2 stats file, returns dict with values.

        Args:
            stats_file (str): Bowtie2 output file with alignment statistics.
        """
        stats = {
            "readCount": None,
            "unpaired": None,
            "unaligned": None,
            "unique": None,
            "multiple": None,
            "alignmentRate": None,
        }
        try:
            with open(stats_file) as handle:
                content = handle.readlines()  # list of strings per line
        except Exception:
            return stats
        # total reads
        try:
            line = [i for i in range(len(content)) if " reads; of these:" in content[i]][0]
            stats["readCount"] = re.sub(r"\D.*", "", content[line])
            if 7 > len(content) > 2:
                line = [
                    i for i in range(len(content)) if "were unpaired; of these:" in content[i]
                ][0]
                stats["unpaired"] = re.sub(r"\D", "", re.sub(r"\(.*", "", content[line]))
            else:
                line = [i for i in range(len(content)) if "were paired; of these:" in content[i]][
                    0
                ]
                stats["unpaired"] = stats["readCount"] - int(
                    re.sub(r"\D", "", re.sub(r"\(.*", "", content[line]))
                )
            line = [i for i in range(len(content)) if "aligned 0 times" in content[i]][0]
            stats["unaligned"] = re.sub(r"\D", "", re.sub(r"\(.*", "", content[line]))
            line = [i for i in range(len(content)) if "aligned exactly 1 time" in content[i]][0]
            stats["unique"] = re.sub(r"\D", "", re.sub(r"\(.*", "", content[line]))
            line = [i for i in range(len(content)) if "aligned >1 times" in content[i]][0]
            stats["multiple"] = re.sub(r"\D", "", re.sub(r"\(.*", "", content[line]))
            line = [i for i in range(len(content)) if "overall alignment rate" in content[i]][0]
            stats["alignmentRate"] = re.sub(r"\%.*", "", content[line]).strip()
        except IndexError:
            pass
        return stats

    def parse_duplicate_stats(self, stats_file: str) -> dict:
        """
        Parses sambamba markdup output, returns dict with values.

        Args:
            stats_file (str): sambamba output file with duplicate statistics.
        """
        series = {}
        try:
            with open(stats_file) as handle:
                content = handle.readlines()  # list of strings per line
        except Exception:
            return series
        try:
            line = [i for i in range(len(content)) if "single ends (among them " in content[i]][0]
            series["single-ends"] = re.sub(r"\D", "", re.sub(r"\(.*", "", content[line]))
            line = [i for i in range(len(content)) if " end pairs...   done in " in content[i]][0]
            series["paired-ends"] = re.sub(r"\D", "", re.sub(r"\.\.\..*", "", content[line]))
            line = [
                i
                for i in range(len(content))
                if " duplicates, sorting the list...   done in " in content[i]
            ][0]
            series["duplicates"] = re.sub(r"\D", "", re.sub(r"\.\.\..*", "", content[line]))
        except IndexError:
            pass
        return series

    def parse_qc(self, qc_file: str) -> dict:
        """
        Parse phantompeakqualtools (spp) QC table and return quality metrics.

        Args:
            qc_file (str): Path to phantompeakqualtools output file, which
                contains sample quality measurements.
        """
        series = {}
        try:
            with open(qc_file) as handle:
                line = handle.readlines()[0].strip().split("\t")  # list of strings per line
            series["NSC"] = line[-3]
            series["RSC"] = line[-2]
            series["qualityTag"] = line[-1]
        except Exception:
            pass
        return series

    def get_peak_number(self, sample: Any) -> Any:
        """
        Counts number of peaks from a sample's peak file.

        Args:
            sample (pipelines.Sample): Sample object with "peaks" attribute.
        """
        proc = subprocess.Popen(["wc", "-l", sample.peaks], stdout=subprocess.PIPE)
        out, err = proc.communicate()
        sample["peakNumber"] = re.sub(r"\D.*", "", out)
        return sample

    def get_frip(self, sample: Any) -> dict:
        """
        Calculates the fraction of reads in peaks for a given sample.

        Args:
            sample (pipelines.Sample): Sample object with "peaks" attribute.
        """
        with open(sample.frip, "r") as handle:
            content = handle.readlines()
        reads_in_peaks = int(re.sub(r"\D", "", content[0]))
        mapped_reads = sample["readCount"] - sample["unaligned"]
        return {"FRiP": reads_in_peaks / mapped_reads}
