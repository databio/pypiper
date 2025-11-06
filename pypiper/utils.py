"""Shared utilities"""

import os
import re
import sys

from collections.abc import Iterable, Mapping, Sequence
from shlex import split
from subprocess import PIPE

if sys.version_info < (3,):
    CHECK_TEXT_TYPES = (str, unicode)
    from inspect import getargspec as get_fun_sig
else:
    CHECK_TEXT_TYPES = (str,)
    from inspect import getfullargspec as get_fun_sig

from ubiquerg import expandpath, is_command_callable

from .const import (
    CHECKPOINT_EXTENSION,
    PIPELINE_CHECKPOINT_DELIMITER,
    STAGE_NAME_SPACE_REPLACEMENT,
)
from .flags import FLAGS

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


# What to export/attach to pypiper package namespace.
# Conceptually, reserve this for functions expected to be used in other
# packages, and import from utils within pypiper for other functions.
__all__ = [
    "add_pypiper_args",
    "build_command",
    "check_all_commands",
    "determine_uncallable",
    "get_first_value",
    "head",
    "logger_via_cli",
    "result_formatter_markdown",
]


CHECKPOINT_SPECIFICATIONS = ["start_point", "stop_before", "stop_after"]


def add_pypiper_args(
    parser, groups=("pypiper",), args=None, required=None, all_args=False
):
    """Use this to add standardized pypiper arguments to your python pipeline.

    There are two ways to use `add_pypiper_args`: by specifying argument groups,
    or by specifying individual arguments. Specifying argument groups will add
    multiple arguments to your parser; these convenient argument groupings
    make it easy to add arguments to certain types of pipeline. For example,
    to make a looper-compatible pipeline, use `groups = ["pypiper", "looper"]`.

    Args:
        parser: ArgumentParser object from a pipeline.
        groups: Adds arguments belonging to specified group
            of args. Options: pypiper, config, looper, resources, common, ngs, all.
        args: You may specify a list of specific arguments one by one.
        required: Arguments to be flagged as 'required' by argparse.
        all_args: Whether to include all of pypiper's arguments defined here.

    Returns:
        A new ArgumentParser object, with selected pypiper arguments added.
    """
    args_to_add = _determine_args(
        argument_groups=groups, arguments=args, use_all_args=all_args
    )
    parser = _add_args(parser, args_to_add, required)
    return parser


def build_command(chunks):
    """Create a command from various parts.

    The parts provided may include a base, flags, option-bound arguments, and
    positional arguments. Each element must be either a string or a two-tuple.
    Raw strings are interpreted as either the command base, a pre-joined
    pair (or multiple pairs) of option and argument, a series of positional
    arguments, or a combination of those elements. The only modification they
    undergo is trimming of any space characters from each end.

    Args:
        chunks: The collection of the command components to interpret, modify,
            and join to create a single meaningful command.

    Returns:
        The single meaningful command built from the given components.

    Raises:
        ValueError: If no command parts are provided.
    """

    if not chunks:
        raise ValueError("No command parts: {} ({})".format(chunks, type(chunks)))

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
    """Ensure existence of folders for a Sample.

    Args:
        sample: Sample (or instance supporting get() that stores folders paths
            in a 'paths' key, in which the value is a mapping from path name to
            actual folder path).
    """
    for path_name, path in sample.paths.items():
        print("{}: '{}'".format(path_name, path))
        base, ext = os.path.splitext(path)
        if ext:
            print("Skipping file-like: '[}'".format(path))
        elif not os.path.isdir(base):
            os.makedirs(base)


def checkpoint_filename(checkpoint, pipeline_name=None):
    """Translate a checkpoint to a filename.

    This not only adds the checkpoint file extension but also standardizes the
    way in which checkpoint names are mapped to filenames.

    Args:
        checkpoint: Name of a pipeline phase/stage.
        pipeline_name: Name of pipeline to prepend to the checkpoint
            filename; this differentiates checkpoint files, e.g. within the
            same sample output folder but associated with different pipelines,
            in case of the (somewhat probable) scenario of a stage name
            collision between pipelines that processed the same sample and
            wrote to the same output folder.

    Returns:
        Standardized checkpoint name for file, plus extension; null if the
        input is a Stage that's designated as a non-checkpoint.
    """
    # Allow Stage as type for checkpoint parameter's argument without
    # needing to import here the Stage type from stage.py module.
    try:
        base = checkpoint.checkpoint_name
    except AttributeError:
        base = translate_stage_name(checkpoint)
    if pipeline_name:
        base = "{}{}{}".format(pipeline_name, PIPELINE_CHECKPOINT_DELIMITER, base)
    return base + CHECKPOINT_EXTENSION


def checkpoint_filepath(checkpoint, pm):
    """Create filepath for indicated checkpoint.

    Args:
        checkpoint: Pipeline phase/stage or one's name.
        pm: Manager of a pipeline instance, relevant for output folder path.

    Returns:
        Standardized checkpoint name for file, plus extension.

    Raises:
        ValueError: If the checkpoint is given as absolute path that does
            not point within pipeline output folder.
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
                    "folder '{}'".format(checkpoint, pm.outfolder)
                )
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


def check_shell(cmd, shell=None):
    """Determine whether a command appears to involve shell process(es).

    The shell argument can be used to override the result of the check.

    Args:
        cmd: Command to investigate.
        shell: Override the result of the check with this value.

    Returns:
        Whether the command appears to involve shell process(es).
    """
    if isinstance(shell, bool):
        return shell
    return "|" in cmd or ">" in cmd or r"*" in cmd


def check_shell_asterisk(cmd):
    """Determine whether a command appears to involve shell stars.

    Args:
        cmd: Command to investigate.

    Returns:
        Whether the command appears to involve shell stars.
    """
    return r"*" in cmd


def check_all_commands(
    cmds,
    get_bad_result=lambda bads: Exception(
        "{} uncallable commands: {}".format(len(bads), bads)
    ),
    handle=None,
):
    """Determine whether all commands are callable.

    Args:
        cmds: Collection of commands to check for callability.
        get_bad_result: How to create result when at least one command is uncallable.
        handle: How to handle with result of failed check.

    Returns:
        Whether all commands given appear callable.

    Raises:
        TypeError: If error handler is provided but isn't callable or
            isn't a single-argument function.
    """
    bads = determine_uncallable(cmds)
    if not bads:
        return True
    if handle is None:

        def handle(res):
            if isinstance(res, Exception):
                raise res
            print("Command check result: {}".format(res))

    elif not hasattr(handle, "__call__") or not 1 == len(get_fun_sig(handle).args):
        raise TypeError("Command check error handler must be a one-arg function")
    handle(get_bad_result(bads))
    return False


def determine_uncallable(
    commands,
    transformations=(
        (
            lambda f: isinstance(f, str)
            and os.path.isfile(expandpath(f))
            and expandpath(f).endswith(".jar"),
            lambda f: "java -jar {}".format(expandpath(f)),
        ),
    ),
    accumulate=False,
):
    """Determine which commands are not callable.

    Args:
        commands: Commands to check for callability.
        transformations: Pairs in which first element is a predicate and second is a
            transformation to apply to the input if the predicate is satisfied.
        accumulate: Whether to accumulate transformations (more than one possible
            per command).

    Returns:
        Collection of commands that appear uncallable; each element is a pair in
        which the first element is the original 'command' and the second is what
        was actually assessed for callability.

    Raises:
        TypeError: If transformations are provided but argument is a string or is
            non-Iterable.
        Exception: If accumulation of transformation is False but the collection of
            transformations is unordered.
    """
    commands = [commands] if isinstance(commands, str) else commands
    if transformations:
        trans = (
            transformations.values()
            if isinstance(transformations, Mapping)
            else transformations
        )
        if (
            not isinstance(transformations, Iterable)
            or isinstance(transformations, str)
            or not all(
                map(
                    lambda func_pair: isinstance(func_pair, tuple)
                    and len(func_pair) == 2,
                    trans,
                )
            )
        ):
            raise TypeError(
                "Transformations argument should be a collection of pairs; got "
                "{} ({})".format(transformations, type(transformations).__name__)
            )
        if accumulate:

            def finalize(cmd):
                for p, t in transformations:
                    if p(cmd):
                        cmd = t(cmd)
                return cmd

        else:
            if not isinstance(transformations, (tuple, list)):
                raise Exception(
                    "If transformations are unordered, non-accumulation of "
                    "effects may lead to nondeterministic behavior."
                )

            def finalize(cmd):
                print("Transformations: {}".format(transformations))
                for p, t in transformations:
                    if p(cmd):
                        return t(cmd)
                return cmd

    else:
        finalize = lambda cmd: cmd
    return [
        (orig, used)
        for orig, used in map(lambda c: (c, finalize(c)), commands)
        if not is_command_callable(used)
    ]


def split_by_pipes_nonnested(cmd):
    """Split the command by shell pipes, but preserve contents in parentheses and braces.

    Args:
        cmd: Command to investigate.

    Returns:
        List of sub commands to be linked.
    """
    # for posterity, this one will do parens only: re.compile(r'(?:[^|(]|\([^)]*\))+')
    # r = re.compile(r'(?:[^|({]|[\({][^)}]*[\)}])+')
    r = re.compile(r"(?:[^|(]|\([^)]*\)+|\{[^}]*\})")
    return r.findall(cmd)


# cmd="a | b | { c | d } ABC | { x  | y } hello '( () e |f )"


def split_by_pipes(cmd):
    """Split the command by shell pipes, but preserve contents in parentheses and braces.

    Also handles nested parens and braces.

    Args:
        cmd: Command to investigate.

    Returns:
        List of sub commands to be linked.
    """

    # Build a simple finite state machine to split on pipes, while
    # handling nested braces or parentheses.
    stack_brace = []
    stack_paren = []
    cmdlist = []
    newcmd = str()
    for char in cmd:
        if char == "{":
            stack_brace.append("{")
        elif char == "}":
            stack_brace.pop()
        elif char == "(":
            stack_paren.append("(")
        elif char == ")":
            stack_paren.pop()

        if len(stack_brace) > 0 or len(stack_paren) > 0:
            # We are inside a parenthetic of some kind; emit character
            # no matter what it is
            newcmd += char
        elif char == "|":
            # if it's a pipe, finish the command and start a new one
            cmdlist.append(newcmd)
            newcmd = str()
            next
        else:
            # otherwise, emit character.
            newcmd += char

    # collect the final command before returning
    cmdlist.append(newcmd)
    return cmdlist


def check_shell_pipes(cmd):
    """Determine whether a command appears to contain shell pipes.

    Args:
        cmd: Command to investigate.

    Returns:
        Whether the command appears to contain shell pipes.
    """
    return "|" in strip_braced_txt(cmd)


def strip_braced_txt(cmd):
    curly_braces = True
    while curly_braces:
        SRE_match_obj = re.search(r"\{(.*?)}", cmd)
        if not SRE_match_obj is None:
            cmd = cmd[: SRE_match_obj.start()] + cmd[(SRE_match_obj.end() + 1) :]
            if re.search(r"\{(.*?)}", cmd) is None:
                curly_braces = False
        else:
            curly_braces = False

    return cmd


def check_shell_redirection(cmd):
    """Determine whether a command appears to contain shell redirection symbol outside of curly brackets.

    Args:
        cmd: Command to investigate.

    Returns:
        Whether the command appears to contain shell redirection.
    """

    return ">" in strip_braced_txt(cmd)


def clear_flags(pm, flag_names=None):
    """Clear pipeline flags.

    Args:
        pm: Pipeline or PipelineManager for which to remove flags.
        flag_names: Names of flags to remove, optional; if unspecified, all
            known flag names will be used.

    Returns:
        Collection of names of flags removed.
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
    """Determine the name for a flag file of the status indicated.

    Args:
        status: Name of status for which to create flag file name.

    Returns:
        Name of flag file corresponding to given status.
    """
    return status + ".flag"


def get_proc_name(cmd):
    """Get the representative process name from complex command.

    Args:
        cmd: A command to be processed.

    Returns:
        The basename representative command.
    """
    if isinstance(cmd, Iterable) and not isinstance(cmd, str):
        cmd = " ".join(cmd)

    return cmd.split()[0].replace("(", "").replace(")", "")


def get_first_value(param, param_pools, on_missing=None, error=True):
    """Get the value for a particular parameter from the first pool in the provided priority list.

    Args:
        param: Name of parameter for which to determine/fetch value.
        param_pools: Ordered (priority) collection of mapping from parameter name
            to value; this should be ordered according to descending priority.
        on_missing: Default value or action to take if the requested parameter is
            missing from all of the pools. If a callable, it should return a value
            when passed the requested parameter as the one and only argument.
        error: Whether to raise an error if the requested parameter is not mapped
            to a value AND there's no value or strategy provided with 'on_missing'
            with which to handle the case of a request for an unmapped parameter.

    Returns:
        Value to which the requested parameter first mapped in the (descending)
        priority collection of parameter 'pools,' or a value explicitly defined
        or derived with 'on_missing.'

    Raises:
        KeyError: If the requested parameter is unmapped in all of the provided
            pools, and the argument to the 'error' parameter evaluates to True.
    """

    # Search for the requested parameter.
    for pool in param_pools:
        if param in pool:
            return pool[param]

    # Raise error if not found and no strategy or value is provided or handling
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
                "a value."
            )
        return on_missing


def head(obj):
    """Facilitate syntactic uniformity for notion of "object-or-first-item."

    Args:
        obj: Single item or collection of items.

    Returns:
        Object itself if non-Iterable, otherwise the first element if one exists.

    Raises:
        ValueError: If the given object is an empty Iterable.
    """
    if isinstance(obj, CHECK_TEXT_TYPES):
        return obj
    try:
        return next(iter(obj))
    except TypeError:
        return obj
    except StopIteration:
        raise ValueError("Requested head of empty iterable")


def is_in_file_tree(fpath, folder):
    """Determine whether a file is in a folder.

    Args:
        fpath: Filepath to investigate.
        folder: Path to folder to query.

    Returns:
        Whether the path indicated is in the folder indicated.
    """
    file_folder, _ = os.path.split(fpath)
    other_folder = os.path.join(folder, "")
    return other_folder.startswith(file_folder)


def is_fastq(file_name):
    """Determine whether indicated file appears to be in FASTQ format.

    Args:
        file_name: Name/path of file to check as FASTQ.

    Returns:
        Whether indicated file appears to be in FASTQ format, zipped or unzipped.
    """
    return is_unzipped_fastq(file_name) or is_gzipped_fastq(file_name)


def is_gzipped_fastq(file_name):
    """Determine whether indicated file appears to be a gzipped FASTQ.

    Args:
        file_name: Name/path of file to check as gzipped FASTQ.

    Returns:
        Whether indicated file appears to be in gzipped FASTQ format.
    """
    _, ext = os.path.splitext(file_name)
    return file_name.endswith(".fastq.gz") or file_name.endswith(".fq.gz")


def is_unzipped_fastq(file_name):
    """Determine whether indicated file appears to be an unzipped FASTQ.

    Args:
        file_name: Name/path of file to check as unzipped FASTQ.

    Returns:
        Whether indicated file appears to be in unzipped FASTQ format.
    """
    _, ext = os.path.splitext(file_name)
    return ext in [".fastq", ".fq"]


def is_sam_or_bam(file_name):
    """Determine whether a file appears to be in a SAM format.

    Args:
        file_name: Name/path of file to check as SAM-formatted.

    Returns:
        Whether file appears to be SAM-formatted.
    """
    _, ext = os.path.splitext(file_name)
    return ext in [".bam", ".sam"]


def logger_via_cli(opts, **kwargs):
    """Build and initialize logger from CLI specification.

    Args:
        opts: Parse of command-line interface.
        **kwargs: Keyword arguments to pass along to underlying logmuse function.

    Returns:
        Newly created and configured logger.
    """
    from copy import deepcopy

    import logmuse

    kwds = deepcopy(kwargs)
    # By default, don't require the logging options to have been added to the parser.
    kwds.setdefault("strict", False)
    return logmuse.logger_via_cli(opts, **kwds)


def make_lock_name(original_path, path_base_folder):
    """Create name for lock file from an absolute path.

    The original path must be absolute, and it should point to a location
    within the location indicated by the base folder path provided. This is
    particularly useful for deleting a sample's output folder path from
    within the path of a target file to generate a lock file corresponding
    to the original target.

    Args:
        original_path: Full original filepath.
        path_base_folder: Portion of original path to delete.

    Returns:
        Name or perhaps relative (to the base folder path indicated) path to lock file.
    """

    def make_name(p):
        if p:
            return p.replace(path_base_folder, "").replace(os.sep, "__")
        else:
            return None

    if isinstance(original_path, str):
        return make_name(original_path)
    elif isinstance(original_path, Sequence):
        result = [make_name(p) for p in original_path]
        return [x for x in result if x]
    raise TypeError(
        "Neither string nor other sequence type: {} ({})".format(
            original_path, type(original_path)
        )
    )


def is_multi_target(target):
    """Determine if pipeline manager's run target is multiple.

    Args:
        target: 0, 1, or multiple targets.

    Returns:
        Whether there are multiple targets.

    Raises:
        TypeError: If the argument is neither None nor string nor Sequence.
    """
    if target is None or isinstance(target, str):
        return False
    elif isinstance(target, Sequence):
        return len(target) > 1
    else:
        raise TypeError(
            "Could not interpret argument as a target: {} ({})".format(
                target, type(target)
            )
        )


def parse_cmd(cmd, shell):
    """Create a list of Popen-distable dicts of commands.

    The commands are split by pipes, if possible.

    Args:
        cmd: The command.
        shell: If the command should be run in the shell rather that in a subprocess.

    Returns:
        List of dicts of commands.
    """

    def _make_dict(command):
        a, s = (
            (command, True) if check_shell(command, shell) else (split(command), False)
        )
        return dict(args=a, stdout=PIPE, shell=s)

    return (
        [_make_dict(c) for c in split_by_pipes(cmd)]
        if not shell and check_shell_pipes(cmd)
        else [dict(args=cmd, stdout=None, shell=True)]
    )


def parse_cores(cores, pm, default):
    """Framework to finalize number of cores for an operation.

    Some calls to a function may directly provide a desired number of cores,
    others may not. Similarly, some pipeline managers may define a cores count
    while others will not. This utility provides a single via which the
    count of cores to use for an operation may be determined. If a cores
    count is given explicitly, use that. Then try pipeline manager for cores.
    Finally, fall back to a default. Force default to be defined (this
    function is intended to be partially applied, then reused within a
    module, class, etc. to standardize the way in which this value is
    determined within a scope.)

    Args:
        cores: Direct specification of cores count.
        pm: Pipeline manager perhaps defining cores.
        default: Default number of cores, used if a value isn't directly given
            and the pipeline manager doesn't define core count.

    Returns:
        Number of cores.
    """
    cores = cores or getattr(pm, "cores", default)
    return int(cores)


def parse_stage_name(stage):
    """Determine the name of a stage.

    The stage may be provided already as a name, as a Stage object, or as a
    callable with __name__ (e.g., function).

    Args:
        stage: Object representing a stage, from which to obtain name.

    Returns:
        Name of putative pipeline Stage.
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
    """Derive path to file for managed pipeline.

    Args:
        pm: Manager of a particular pipeline instance.
        filename: Name of file for which to create full path based on pipeline's
            output folder.
        suffix: Suffix for the file; this can be added to the filename if provided
            or added to the pipeline name if there's no filename.

    Returns:
        Path to file within managed pipeline's output folder, with filename as
        given or determined by the pipeline name, and suffix appended if given.

    Raises:
        TypeError: If neither filename nor suffix is provided, raise a TypeError,
            as in that case there's no substance from which to create a filepath.
    """
    if filename is None and suffix is None:
        raise TypeError(
            "Provide filename and/or suffix to create path to a pipeline file."
        )
    filename = (filename or pm.name) + (suffix or "")

    # Note that Pipeline and PipelineManager define the same outfolder.
    # In fact, a Pipeline just references its manager's outfolder.
    # So we can handle argument of either type to pm parameter.
    return filename if os.path.isabs(filename) else os.path.join(pm.outfolder, filename)


def translate_stage_name(stage):
    """Account for potential variability in stage/phase name definition.

    Since a pipeline author is free to name his/her processing phases/stages
    as desired, but these choices influence file names, enforce some
    standardization. Specifically, prohibit potentially problematic spaces.

    Args:
        stage: Pipeline stage, its name, or a representative function.

    Returns:
        Standardized pipeline phase/stage name.
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
    """Determine the arguments to add to a parser (for a pipeline).

    Args:
        argument_groups: Collection of names of groups of arguments to add to an
            argument parser.
        arguments: Collection of specific arguments to add to the parser.
        use_all_args: Whether to use all arguments defined here.

    Returns:
        Collection of (unique) argument names to add to a parser.
    """

    from collections.abc import Iterable

    from logmuse import LOGGING_CLI_OPTDATA

    # Define the argument groups.
    args_by_group = {
        "pypiper": ["recover", "new-start", "dirty", "force-follow", "testmode"]
        + [*LOGGING_CLI_OPTDATA],
        "config": ["config"],
        "checkpoint": ["stop-before", "stop-after"],
        "resource": ["mem", "cores"],
        "looper": ["config", "output-parent", "mem", "cores", "pipeline-name"],
        "common": ["input", "sample-name"],
        "ngs": ["sample-name", "input", "input2", "genome", "single-or-paired"],
        "logmuse": [*LOGGING_CLI_OPTDATA],
        "pipestat": [
            "pipestat-namespace",
            "pipestat-record-id",
            "pipestat-schema",
            "pipestat-results-file",
            "pipestat-config",
        ],
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


def default_pipeline_config(pipeline_filepath):
    """Determine the default filepath for a pipeline's config file.

    Args:
        pipeline_filepath: Path to a pipeline.

    Returns:
        Default filepath for pipeline's config file.
    """
    return os.path.splitext(os.path.basename(pipeline_filepath))[0] + ".yaml"


def default_pipestat_output_schema(pipeline_filepath):
    """Determine the default filepath for a pipeline's pipestat output schema.

    Args:
        pipeline_filepath: Path to a pipeline.

    Returns:
        Default filepath for a pipeline's pipestat output schema.
    """
    pipestat_results_schema = os.path.join(
        os.path.dirname(pipeline_filepath), "pipestat_output_schema.yaml"
    )
    print(f"Using default schema: {pipestat_results_schema}")
    return pipestat_results_schema if os.path.exists(pipestat_results_schema) else None


def _add_args(parser, args, required):
    """Add new arguments to an ArgumentParser.

    Args:
        parser: Instance to update with new arguments.
        args: Collection of names of arguments to add.
        required: Collection of arguments to designate as required.

    Returns:
        Updated ArgumentParser.
    """

    import copy

    required = required or []

    default_config = default_pipeline_config(sys.argv[0])

    # Define the arguments.
    argument_data = {
        "testmode": (
            "-T",
            {"action": "store_true", "help": "Only print commands, don't run"},
        ),
        "recover": (
            "-R",
            {
                "action": "store_true",
                "help": "Overwrite locks to recover from previous failed run",
            },
        ),
        "new-start": (
            "-N",
            {
                "action": "store_true",
                "help": "Overwrite all results to start a fresh run",
            },
        ),
        "dirty": (
            "-D",
            {"action": "store_true", "help": "Don't auto-delete intermediate files"},
        ),
        "force-follow": (
            "-F",
            {"action": "store_true", "help": "Always run 'follow' commands"},
        ),
        "start-point": {"help": "Name of pipeline stage at which to begin"},
        "stop-before": {
            "help": "Name of pipeline stage at which to stop "
            "(exclusive, i.e. not run)"
        },
        "stop-after": {
            "help": "Name of pipeline stage at which to stop " "(inclusive, i.e. run)"
        },
        "config": (
            "-C",
            {
                "dest": "config_file",
                "metavar": "CONFIG_FILE",
                "default": default_config,
                "help": "Pipeline configuration file (YAML). "
                "Relative paths are with respect to the "
                "pipeline script.",
            },
        ),
        "pipeline-name": {"metavar": "PIPELINE_NAME", "help": "Name of the pipeline"},
        "sample-name": (
            "-S",
            {"metavar": "SAMPLE_NAME", "help": "Name for sample to run"},
        ),
        "output-parent": (
            "-O",
            {
                "metavar": "PARENT_OUTPUT_FOLDER",
                "help": "Parent output directory of project",
            },
        ),
        "cores": (
            "-P",
            {
                "type": int,
                "default": 1,
                "metavar": "NUMBER_OF_CORES",
                "help": "Number of cores for parallelized processes",
            },
        ),
        "mem": (
            "-M",
            {
                "default": "4000",
                "metavar": "MEMORY_LIMIT",
                "help": "Memory limit for processes accepting such. "
                "Default units are megabytes unless specified "
                "using the suffix [K|M|G|T].",
            },
        ),
        "input": (
            "-I",
            {
                "nargs": "+",
                "metavar": "INPUT_FILES",
                "help": "One or more primary input files",
            },
        ),
        "input2": (
            "-I2",
            {
                "nargs": "*",
                "metavar": "INPUT_FILES2",
                "help": "Secondary input files, such as read2",
            },
        ),
        "genome": (
            "-G",
            {"dest": "genome_assembly", "help": "Identifier for genome assembly"},
        ),
        "single-or-paired": (
            "-Q",
            {"default": "single", "help": "Single- or paired-end sequencing protocol"},
        ),
        "pipestat-namespace": {
            "help": "Namespace to report into. This will be the DB table name "
            "if using DB as the object back-end"
        },
        "pipestat-record-id": {"help": "Record identifier to report for"},
        "pipestat-schema": {
            "help": "Path to the output schema that formalizes the " "results structure"
        },
        "pipestat-config": {"help": "Path to the configuration file"},
        "pipestat-results-file": {
            "help": "YAML file to report into, if file is used as "
            "the object back-end"
        },
    }

    from logmuse import LOGGING_CLI_OPTDATA

    argument_data.update(LOGGING_CLI_OPTDATA)

    if len(required) > 0:
        required_named = parser.add_argument_group("required named arguments")

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
                    "specification data."
                )

        argdata["required"] = arg in required

        long_opt = "--{}".format(arg)
        opts = (short_opt, long_opt) if short_opt else (long_opt,)
        if arg in required:
            required_named.add_argument(*opts, **argdata)
        else:
            parser.add_argument(*opts, **argdata)

    return parser


def result_formatter_markdown(pipeline_name, record_identifier, res_id, value) -> str:
    """Returns Markdown formatted value as string.

    Note: Pipeline_name and record_identifier should be kept because pipestat needs it.

    Args:
        pipeline_name: Name of the pipeline.
        record_identifier: Identifier for the record.
        res_id: Result identifier.
        value: Value to format.

    Returns:
        Markdown formatted value as string.
    """

    message_markdown = "\n> `{key}`\t{value}\t_RES_".format(key=res_id, value=value)

    return message_markdown
