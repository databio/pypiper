""" Shared utilities """

from collections import Iterable
import os
import sys
import re

from .const import \
    CHECKPOINT_EXTENSION, PIPELINE_CHECKPOINT_DELIMITER, \
    STAGE_NAME_SPACE_REPLACEMENT
from .flags import FLAGS


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


# What to export/attach to pypiper package namespace.
# Conceptually, reserve this for functions expected to be used in other
# packages, and import from utils within pypiper for other functions.
__all__ = ["add_pypiper_args", "build_command", "get_first_value"]


CHECKPOINT_SPECIFICATIONS = ["start_point", "stop_before", "stop_after"]


def add_pypiper_args(parser, groups=("pypiper", ), args=None,
                     required=None, all_args=False):
    """
    Use this to add standardized pypiper arguments to your python pipeline.

    There are two ways to use `add_pypiper_args`: by specifying argument groups,
    or by specifying individual arguments. Specifying argument groups will add
    multiple arguments to your parser; these convenient argument groupings
    make it easy to add arguments to certain types of pipeline. For example,
    to make a looper-compatible pipeline, use `groups = ["pypiper", "looper"]`.

    :param parser: an ArgumentParser object from your pipeline
    :type parser: argparse.ArgumentParser
    :param groups: Adds arguments belong to specified group of args.
         Options are: pypiper, config, looper, resources, common, ngs, all.
    :type groups: Iterable[str] | str
    :param args: You may specify a list of specific arguments one by one.
    :type args: Iterable[str] | str
    :param required: Arguments to be flagged as 'required' by argparse.
    :type required: Iterable[str]
    :param all_args: Whether to include all of pypiper's arguments defined here.
    :type all_args: bool
    :return: A new ArgumentParser object, with selected pypiper arguments added
    :rtype: argparse.ArgumentParser
    """
    args_to_add = _determine_args(
        argument_groups=groups, arguments=args, use_all_args=all_args)
    parser = _add_args(parser, args_to_add, required)
    return parser


def build_command(chunks):
    """
    Create a command from various parts.

    The parts provided may include a base, flags, option-bound arguments, and
    positional arguments. Each element must be either a string or a two-tuple.
    Raw strings are interpreted as either the command base, a pre-joined
    pair (or multiple pairs) of option and argument, a series of positional
    arguments, or a combination of those elements. The only modification they
    undergo is trimming of any space characters from each end.

    :param Iterable[str | (str, str | NoneType)] chunks: the collection of the
        command components to interpret, modify, and join to create a
        single meaningful command
    :return str: the single meaningful command built from the given components
    :raise ValueError: if no command parts are provided
    """

    if not chunks:
        raise ValueError(
            "No command parts: {} ({})".format(chunks, type(chunks)))

    if isinstance(chunks, str):
        return chunks

    parsed_pieces = []

    for cmd_part in chunks:
        if cmd_part is None:
            continue
        try:
            # Trim just space, not all whitespace.
            # This prevents damage to an option that specifies,
            # say, tab as a delimiter.
            parsed_pieces.append(cmd_part.strip(" "))
        except AttributeError:
            option, argument = cmd_part
            if argument is None or argument == "":
                continue
            option, argument = option.strip(" "), str(argument).strip(" ")
            parsed_pieces.append("{} {}".format(option, argument))

    return " ".join(parsed_pieces)


def build_sample_paths(sample):
    """
    Ensure existence of folders for a Sample.

    :param looper.models.Sample sample: Sample (or instance supporting get()
        that stores folders paths in a 'paths' key, in which the value is a
        mapping from path name to actual folder path)
    """
    for path_name, path in sample.paths.items():
        print("{}: '{}'".format(path_name, path))
        base, ext = os.path.splitext(path)
        if ext:
            print("Skipping file-like: '[}'".format(path))
        elif not os.path.isdir(base):
            os.makedirs(base)


def checkpoint_filename(checkpoint, pipeline_name=None):
    """
    Translate a checkpoint to a filename.

    This not only adds the checkpoint file extension but also standardizes the
    way in which checkpoint names are mapped to filenames.

    :param checkpoint: name of a pipeline phase/stage
    :type checkpoint: str | Stage
    :param pipeline_name: name of pipeline to prepend to the checkpoint
        filename; this differentiates checkpoint files, e.g. within the
        same sample output folder but associated with different pipelines,
        in case of the (somewhat probable) scenario of a stage name
        collision between pipelines the processed the same sample and
        wrote to the same output folder
    :type pipeline_name: str
    :return: standardized checkpoint name for file, plus extension;
        null if the input is a Stage that's designated as a non-checkpoint
    :rtype: str | NoneType
    """
    # Allow Stage as type for checkpoint parameter's argument without
    # needing to import here the Stage type from stage.py module.
    try:
        base = checkpoint.checkpoint_name
    except AttributeError:
        base = translate_stage_name(checkpoint)
    if pipeline_name:
        base = "{}{}{}".format(
            pipeline_name, PIPELINE_CHECKPOINT_DELIMITER, base)
    return base + CHECKPOINT_EXTENSION


def checkpoint_filepath(checkpoint, pm):
    """
    Create filepath for indicated checkpoint.

    :param checkpoint: Pipeline phase/stage or one's name
    :type checkpoint: str | Stage
    :param pm: manager of a pipeline instance, relevant for output folder path.
    :type pm: pypiper.PipelineManager | pypiper.Pipeline
    :return str: standardized checkpoint name for file, plus extension
    :raise ValueError: if the checkpoint is given as absolute path that does
        not point within pipeline output folder
    """

    # Handle case in which checkpoint is given not just as a string, but
    # as a checkpoint-like filename. Don't worry about absolute path status
    # of a potential filename input, or whether it's in the pipeline's
    # output folder. That's handled upstream. While this isn't a protected
    # function, there's no real reason to call this from outside the package.
    if isinstance(checkpoint, str):
        if os.path.isabs(checkpoint):
            if is_in_file_tree(checkpoint, pm.outfolder):
                return checkpoint
            else:
                raise ValueError(
                    "Absolute checkpoint path '{}' is not in pipeline output "
                    "folder '{}'".format(checkpoint, pm.outfolder))
        _, ext = os.path.splitext(checkpoint)
        if ext == CHECKPOINT_EXTENSION:
            return pipeline_filepath(pm, filename=checkpoint)

    # Allow Pipeline as pm type without importing Pipeline.
    try:
        pm = pm.manager
    except AttributeError:
        pass

    # We want the checkpoint filename itself to become a suffix, with a
    # delimiter intervening between the pipeline name and the checkpoint
    # name + extension. This is to handle the case in which a single, e.g.,
    # sample's output folder is the destination for output from multiple
    # pipelines, and we thus want to be able to distinguish between
    # checkpoint files from different pipelines for that sample that may
    # well define one or more stages with the same name (e.g., trim_reads,
    # align_reads, etc.)
    chkpt_name = checkpoint_filename(checkpoint, pipeline_name=pm.name)
    return pipeline_filepath(pm, filename=chkpt_name)


def check_shell(cmd):
    """
    Determine whether a command appears to involve shell process(es).

    :param str cmd: Command to investigate.
    :return bool: Whether the command appears to involve shell process(es).
    """
    return "|" in cmd or ">" in cmd or r"*" in cmd


def check_shell_asterisk(cmd):
    """
    Determine whether a command appears to involve shell stars.

    :param str cmd: Command to investigate.
    :return bool: Whether the command appears to involve shell stars.
    """
    return r"*" in cmd

def split_by_pipes(cmd):
    """
    Split the command by shell pipes, but preserve the contents of the parantheses.

    :param str cmd: Command to investigate.
    :return list: List of sub commands to be linked
    """
    r = re.compile(r'(?:[^|(]|\([^)]*\))+')
    return r.findall(cmd)

def check_shell_pipes(cmd):
    """
    Determine whether a command appears to contain shell pipes.

    :param str cmd: Command to investigate.
    :return bool: Whether the command appears to contain shell pipes.
    """
    return "|" in cmd


def check_shell_redirection(cmd):
    """
    Determine whether a command appears to contain shell redirection symbol outside of curly brackets

    :param str cmd: Command to investigate.
    :return bool: Whether the command appears to contain shell redirection.
    """
    curly_brackets = True
    while curly_brackets:
        SRE_match_obj = re.search(r'\{(.*?)}',cmd)
        if SRE_match_obj is not None:
            cmd = cmd[:SRE_match_obj.start()] + cmd[(SRE_match_obj.end()+1):]
            if re.search(r'\{(.*?)}',cmd) is None:
                curly_brackets = False
        else:
            curly_brackets = False
    return ">" in cmd


def clear_flags(pm, flag_names=None):
    """

    :param pm: Pipeline or PipelineManager for which to remove flags
    :type pm: pypiper.PipelineManager | pypiper.Pipeline
    :param flag_names: Names of flags to remove, optional; if unspecified,
        all known flag names will be used.
    :type flag_names: Iterable[str]
    :return: Collection of names of flags removed
    """

    # Accept Pipeline as type of pm without needing to import Pipeline.
    try:
        pm = pm.manager
    except AttributeError:
        pass
    flag_names = flag_names or FLAGS
    if isinstance(flag_names, str):
        flag_names = [flag_names]
    removed = []
    for f in flag_names:
        flag_file_suffix = "_{}".format(flag_name(f))
        path_flag_file = pipeline_filepath(pm, suffix=flag_file_suffix)
        try:
            os.remove(path_flag_file)
        except:
            pass
        else:
            print("Removed existing flag: '{}'".format(path_flag_file))
            removed.append(f)
    return removed


def flag_name(status):
    """
    Determine the name for a flag file of the status indicated.

    :param status: Name of status for which to create flag file name.
    :type status: str
    :return: Name of flag file corresponding to given status.
    :rtype: str
    """
    return status + ".flag"


def get_first_value(param, param_pools, on_missing=None, error=True):
    """
    Get the value for a particular parameter from the first pool in the provided
    priority list of parameter pools.

    :param str param: Name of parameter for which to determine/fetch value.
    :param Sequence[Mapping[str, object]] param_pools: Ordered (priority)
        collection of mapping from parameter name to value; this should be
        ordered according to descending priority.
    :param object | function(str) -> object on_missing: default value or
        action to take if the requested parameter is missing from all of the
        pools. If a callable, it should return a value when passed the
        requested parameter as the one and only argument.
    :param bool error: Whether to raise an error if the requested parameter
        is not mapped to a value AND there's no value or strategy provided
        with 'on_missing' with which to handle the case of a request for an
        unmapped parameter.
    :return object: Value to which the requested parameter first mapped in
        the (descending) priority collection of parameter 'pools,' or
        a value explicitly defined or derived with 'on_missing.'
    :raise KeyError: If the requested parameter is unmapped in all of the
        provided pools, and the argument to the 'error' parameter evaluates
        to True.
    """

    # Search for the requested parameter.
    for pool in param_pools:
        if param in pool:
            return pool[param]

    # Raise error if unfound and no strategy or value is provided or handling
    # unmapped parameter requests.
    if error and on_missing is None:
        raise KeyError("Unmapped parameter: '{}'".format(param))

    # Use the value or strategy for handling unmapped parameter case.
    try:
        return on_missing(param)
    except TypeError:
        if hasattr(on_missing, "__call__"):
            raise TypeError(
                "Any callable passed as the action to take when a requested "
                "parameter is missing should accept that parameter and return "
                "a value.")
        return on_missing


def is_in_file_tree(fpath, folder):
    """
    Determine whether a file is in a folder.

    :param str fpath: filepath to investigate
    :param folder: path to folder to query
    :return bool: whether the path indicated is in the folder indicated
    """
    file_folder, _ = os.path.split(fpath)
    other_folder = os.path.join(folder, "")
    return other_folder.startswith(file_folder)


def is_fastq(file_name):
    """
    Determine whether indicated file appears to be in FASTQ format.

    :param file_name: Name/path of file to check as FASTQ.
    :type file_name: str
    :return: Whether indicated file appears to be in FASTQ format, zipped
        or unzipped.
    :rtype bool
    """
    return is_unzipped_fastq(file_name) or is_gzipped_fastq(file_name)


def is_gzipped_fastq(file_name):
    """
    Determine whether indicated file appears to be a gzipped FASTQ.

    :param file_name: Name/path of file to check as gzipped FASTQ.
    :type file_name: str
    :return: Whether indicated file appears to be in gzipped FASTQ format.
    :rtype bool
    """
    _, ext = os.path.splitext(file_name)
    return file_name.endswith(".fastq.gz") or file_name.endswith(".fq.gz")


def is_unzipped_fastq(file_name):
    """
    Determine whether indicated file appears to be an unzipped FASTQ.

    :param file_name: Name/path of file to check as unzipped FASTQ.
    :type file_name: str
    :return: Whether indicated file appears to be in unzipped FASTQ format.
    :rtype bool
    """
    _, ext = os.path.splitext(file_name)
    return ext in [".fastq", ".fq"]


def is_sam_or_bam(file_name):
    """
    Determine whether a file appears to be in a SAM format.

    :param file_name: Name/path of file to check as SAM-formatted.
    :type file_name: str
    :return: Whether file appears to be SAM-formatted
    :rtype: bool
            """
    _, ext = os.path.splitext(file_name)
    return ext in [".bam", ".sam"]


def make_lock_name(original_path, path_base_folder):
    """
    Create name for lock file from an absolute path.

    The original path must be absolute, and it should point to a location
    within the location indicated by the base folder path provided. This is
    particularly useful for deleting a sample's output folder path from
    within the path of a target file to generate a lock file corresponding
    to the original target.

    :param str original_path: Full original filepath.
    :param str path_base_folder: Portion of original path to delete
    :return str: Name or perhaps relative (to the base folder path indicated)
        path to lock file
    """
    return original_path.replace(path_base_folder, "").replace(os.sep, "__")


def parse_cores(cores, pm, default):
    """
    Framework to finalize number of cores for an operation.

    Some calls to a function may directly provide a desired number of cores,
    others may not. Similarly, some pipeline managers may define a cores count
    while others will not. This utility provides a single via which the
    count of cores to use for an operation may be determined. If a cores
    count is given explicitly, use that. Then try pipeline manager for cores.
    Finally, fall back to a default. Force default to be defined (this
    function is intended to be partially applied, then reused within a
    module, class, etc. to standardize the way in which this value is
    determined within a scope.)

    :param int | str cores: direct specification of cores count
    :param pypiper.PipelineManager pm: pipeline manager perhaps defining cores
    :param int | str default: default number of cores, used if a value isn't
        directly given and the pipeline manager doesn't define core count.
    :return int: number of cores
    """
    cores = cores or getattr(pm, "cores", default)
    return int(cores)


def parse_stage_name(stage):
    """
    Determine the name of a stage.

    The stage may be provided already as a name, as a Stage object, or as a
    callable with __name__ (e.g., function).

    :param stage: Object representing a stage, from which to obtain name.
    :type stage: str | pypiper.Stage | function
    :return: Name of putative pipeline Stage.
    :rtype: str
    """
    if isinstance(stage, str):
        return stage
    try:
        return stage.name
    except AttributeError:
        try:
            return stage.__name__
        except AttributeError:
            raise TypeError("Unsupported stage type: {}".format(type(stage)))


def pipeline_filepath(pm, filename=None, suffix=None):
    """
    Derive path to file for managed pipeline.

    :param pm: Manager of a particular pipeline instance.
    :type pm: pypiper.PipelineManager | pypiper.Pipeline
    :param filename: Name of file for which to create full path based
        on pipeline's output folder.
    :type filename: str
    :param suffix: Suffix for the file; this can be added to the filename
        if provided or added to the pipeline name if there's no filename.
    :type suffix: str
    :raises TypeError: If neither filename nor suffix is provided, raise a
        TypeError, as in that case there's no substance from which to create
        a filepath.
    :return: Path to file within managed pipeline's output folder, with
        filename as given or determined by the pipeline name, and suffix
        appended if given.
    :rtype: str
    """

    if filename is None and suffix is None:
        raise TypeError("Provide filename and/or suffix to create "
                        "path to a pipeline file.")

    filename = (filename or pm.name) + (suffix or "")

    # Note that Pipeline and PipelineManager define the same outfolder.
    # In fact, a Pipeline just references its manager's outfolder.
    # So we can handle argument of either type to pm parameter.
    return filename if os.path.isabs(filename) \
        else os.path.join(pm.outfolder, filename)


def translate_stage_name(stage):
    """
    Account for potential variability in stage/phase name definition.

    Since a pipeline author is free to name his/her processing phases/stages
    as desired, but these choices influence file names, enforce some
    standardization. Specifically, prohibit potentially problematic spaces.

    :param stage: Pipeline stage, its name, or a representative function.
    :type stage: str | pypiper.Stage | function
    :return: Standardized pipeline phase/stage name.
    :rtype: str
    """
    # First ensure that we have text.
    name = parse_stage_name(stage)
    # Cast to string to ensure that indexed stages (ints are handled).
    return str(name).lower().replace(" ", STAGE_NAME_SPACE_REPLACEMENT)


# TODO: implement as context manager.
class Tee(object):
    def __init__(self, log_file):
        self.file = open(log_file, "a")
        self.stdout = sys.stdout
        sys.stdout = self

    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()

    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)
        self.file.flush()
        self.stdout.flush()

    def fileno(self):
        return self.stdout.fileno()


def uniqify(seq):  # Dave Kirby
    # Order preserving
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]


def _determine_args(argument_groups, arguments, use_all_args=False):
    """
    Determine the arguments to add to a parser (for a pipeline).

    :param argument_groups: Collection of names of groups of arguments to
        add to an argument parser.
    :type argument_groups: Iterable[str] | str
    :param arguments: Collection of specific arguments to add to the parser.
    :type arguments: Iterable[str] | str
    :param use_all_args: Whether to use all arguments defined here.
    :type use_all_args: bool
    :return: Collection of (unique) argument names to add to a parser.
    :rtype: Set[str]
    """

    if sys.version_info < (3, 3):
        from collections import Iterable
    else:
        from collections.abc import Iterable

    # Define the argument groups.
    args_by_group = {
        "pypiper": ["recover", "new-start", "dirty", "force-follow"],
        "config": ["config"],
        "checkpoint": ["stop-before", "stop-after"],
        "resource": ["mem", "cores"],
        "looper": ["config", "output-parent", "mem", "cores"],
        "common": ["input", "sample-name"],
        "ngs": ["sample-name", "input", "input2", "genome", "single-or-paired"]
    }

    # Handle various types of group specifications.
    groups = None
    if use_all_args:
        groups = args_by_group.keys()
    elif isinstance(argument_groups, str):
        groups = [argument_groups]
    elif isinstance(argument_groups, Iterable):
        groups = argument_groups
    elif argument_groups:
        raise TypeError("arguments must be a str or a list.")

    # Collect the groups of arguments.
    final_args = list()
    if groups:
        for g in groups:
            try:
                this_group_args = args_by_group[g]
            except KeyError:
                print("Skipping undefined pypiper argument group '{}'".format(g))
            else:
                final_args.extend(this_group_args)
                # final_args |= {this_group_args} if \
                #     isinstance(this_group_args, str) else set(this_group_args)

    # Handle various types of specific, individual argument specifications.
    if isinstance(arguments, str):
        final_args.append(arguments)
    elif isinstance(arguments, Iterable):
        final_args.extend(arguments)
    elif arguments:
        raise TypeError("arguments must be a str or a list.")

    return uniqify(final_args)


def _add_args(parser, args, required):
    """
    Add new arguments to an ArgumentParser.

    :param parser: ArgumentParser to update with new arguments
    :type parser: argparse.ArgumentParser
    :param args: Collection of names of arguments to add.
    :type args: Iterable[str]
    :param required: Collection of arguments to designate as required
    :type required: Iterable[str]
    :return: Updated ArgumentParser
    :rtype: argparse.ArgumentParser
    """

    import copy

    required = required or []

    # Determine the default pipeline config file.
    pipeline_script = os.path.basename(sys.argv[0])
    default_config, _ = os.path.splitext(pipeline_script)
    default_config += ".yaml"

    # Define the arguments.
    argument_data = {
        "recover":
            ("-R", {"action": "store_true",
                    "help": "Overwrite locks to recover from previous failed run"}),
        "new-start":
            ("-N", {"action": "store_true",
                    "help": "Overwrite all results to start a fresh run"}),
        "dirty":
            ("-D", {"action": "store_true",
                    "help": "Don't auto-delete intermediate files"}),
        "force-follow":
            ("-F", {"action": "store_true",
                    "help": "Always run 'follow' commands"}),
        "start-point":
            {"help": "Name of pipeline stage at which to begin"},
        "stop-before":
            {"help": "Name of pipeline stage at which to stop "
                     "(exclusive, i.e. not run)"},
        "stop-after":
            {"help": "Name of pipeline stage at which to stop "
                     "(inclusive, i.e. run)"},
        "config":
            ("-C", {"dest": "config_file", "metavar": "CONFIG_FILE",
                    "default": default_config,
                    "help": "Pipeline configuration file (YAML). "
                            "Relative paths are with respect to the "
                            "pipeline script."}),
        "sample-name":
            ("-S", {"metavar": "SAMPLE_NAME",
                    "help": "Name for sample to run"}),
        "output-parent":
            ("-O", {"metavar": "PARENT_OUTPUT_FOLDER",
                    "help": "Parent output directory of project"}),
        "cores":
            ("-P", {"type": int, "default": 1, "metavar": "NUMBER_OF_CORES",
                    "help": "Number of cores for parallelized processes"}),
        "mem":
            ("-M", {"default": "4000", "metavar": "MEMORY_LIMIT",
                    "help": "Memory limit (in Mb) for processes accepting such"}),
        "input":
            ("-I", {"nargs": "+", "metavar": "INPUT_FILES",
                    "help": "One or more primary input files"}),
        "input2":
            ("-I2", {"nargs": "*", "metavar": "INPUT_FILES2",
                     "help": "Secondary input files, such as read2"}),
        "genome":
            ("-G", {"dest": "genome_assembly",
                    "help": "Identifier for genome assembly"}),
        "single-or-paired":
            ("-Q", {"default": "single",
                    "help": "Single- or paired-end sequencing protocol"})
    }

    if len(required) > 0:
        required_named = parser.add_argument_group('required named arguments')

    # Configure the parser for each argument.
    for arg in args:
        try:
            argdata = copy.deepcopy(argument_data[arg])
        except KeyError:
            print("Skipping undefined pypiper argument: '{}'".format(arg))
            continue
        if isinstance(argdata, dict):
            short_opt = None
        else:
            try:
                short_opt, argdata = argdata
            except ValueError:
                raise TypeError(
                    "Option name must map to dict or two-tuple (short "
                    "name and dict) of argument command-line argument "
                    "specification data.")

        argdata["required"] = arg in required

        long_opt = "--{}".format(arg)
        opts = (short_opt, long_opt) if short_opt else (long_opt, )
        if arg in required:
            required_named.add_argument(*opts, **argdata)
        else:
            parser.add_argument(*opts, **argdata)

    return parser
