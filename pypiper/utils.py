""" Shared utilities """

import os
import sys


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


# What to export/attach to pypiper package namespace.
# Conceptually, reserve this for functions expected to be used in other
# packages, and import from utils within pypiper for other functions.
__all__ = ["add_pypiper_args"]



def add_pypiper_args(parser, groups=("pypiper", ), args=None,
                     required=None, all_args=False):
    """
    Adds default automatic args to an ArgumentParser. Use this to add standardized
    pypiper arguments to your python pipeline.

    There are two ways to use `add_pypiper_args`: by specifying argument groups,
    or by specifying individual arguments. Specifying argument groups will add
    multiple arguments to your parser; these convenient argumenet groupings make it
    easy to add arguments to certain types of pipeline. For example, to make a
    looper-compatible pipeline, use `groups = ["pypiper", "looper"]`.

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



def check_shell(cmd):
    """
    Determine whether a command appears to involve shell proccess(es).

    :param str cmd: Command to investigate.
    :return bool: Whether the command appears to involve shell proccess(es).
    """
    return "|" in cmd or ">" in cmd or r"*" in cmd



def pipeline_filepath(pm, filename=None, suffix=None):
    """
    Derive path to file for managed pipeline.

    :param pm: Manager of a particular pipeline instance.
    :type pm: pypiper.PipelineManager
    :param filename: Name of file for which to create full path based
        on pipeline's output folder.
    :type filename: str
    :param suffix: Suffix for the file; this can be added to the filename
        if provided or added to the pipeline name if there's no filename.
    :type suffix: str
    :return: Path to file within managed pipeline's output folder, with
        filename as given or determined by the pipeline name, and suffix
        appended if given.
    :rtype: str
    :raises ValueError: If neither filename nor suffix is provided, raise a
        ValueError, as in that case there's no substance from which to create
        a filepath.
    """
    if filename is None and suffix is None:
        raise ValueError("Provide filename and/or suffix to create "
                         "path to a pipeline file.")
    filename = (filename or pm.pipeline_name) + (suffix or "")
    return os.path.join(pm.pipeline_outfolder, filename)



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

    # Define the argument groups.
    args_by_group = {
        "pypiper" : ["recover", "new-start", "dirty", "follow"],
        "config" : ["config"],
        "resource" : ["mem", "cores"],
        "looper" : ["config", "output-parent", "mem", "cores"],
        "common" : ["input", "sample-name"],
        "ngs" : ["input", "sample-name", "input2", "genome", "single-or-paired"]
    }

    # Handle various types of group specifications.
    if use_all_args:
        groups = args_by_group.keys()
    elif isinstance(argument_groups, str):
        groups = {argument_groups}
    else:
        groups = set(argument_groups or [])

    # Collect the groups of arguments.
    final_args = set()
    for g in groups:
        try:
            this_group_args = args_by_group[g]
        except KeyError:
            print("Skipping undefined pypiper argument group '{}'".format(g))
        else:
            final_args |= {this_group_args} if \
                isinstance(this_group_args, str) else set(this_group_args)

    # Handle various types of specific, individual argument specifications.
    if isinstance(arguments, str):
        arguments = {arguments}
    else:
        arguments = set(arguments or [])

    return final_args | arguments



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

    # Determine the default pipeline config file.
    pipeline_script = os.path.basename(sys.argv[0])
    default_config, _ = os.path.splitext(pipeline_script)
    default_config += ".yaml"

    # Define the arguments.
    argument_data = {
        "recover":
            ("-R", {"action": "store_true",
                    "help": "Recover mode, overwrite locks"}),
        "new-start":
            ("-N", {"dest": "fresh", "action": "store_true",
                    "help": "Fresh start mode, overwrite all"}),
        "dirty":
            ("-D", {"dest": "manual_clean", "action": "store_true",
                    "help": "Make all cleanups manual"}),
        "follow":
            ("-F", {"dest": "force_follow", "action": "store_true",
                    "help": "Run all 'follow' commands, even if the "
                            "primary command is not run"}),
        "config":
            ("-C", {"dest": "config_file", "metavar": "CONFIG_FILE",
                    "default": default_config,
                    "help": "Pipeline configuration file (YAML). "
                            "Relative paths are with respect to the "
                            "pipeline script."}),
        "output-parent":
            ("-O", {"metavar": "PARENT_OUTPUT_FOLDER",
                    "help": "Parent output directory of project"}),
        "cores":
            ("-P", {"type": int, "default": 1, "metavar": "NUMBER_OF_CORES",
                    "help": "Number of cores for parallelized processes"}),
        "mem":
            ("-M", {"default": "4000", "metavar": "MEMORY_LIMIT",
                    "help": "Amount of memory (Mb) use to allow for "
                            "processes for which that can be specified"}),
        "input":
            ("-I", {"nargs": "+", "metavar": "INPUT_FILES",
                    "help": "One or more primary input files (required)"}),
        "input2":
            ("-I2", {"nargs": "*", "metavar": "INPUT_FILES2",
                     "help": "Secondary input file(s), e.g. read2 for a "
                             "paired-end protocol"}),
        "genome":
            ("-G", {"dest": "genome_assembly",
                    "help": "Identifier for genome assembly"}),
        "single-or-paired":
            ("-Q", {"default": "single",
                    "help": "Single- or paired-end sequencing protocol"})
    }

    # Configure the parser for each argument.
    for arg in args:
        try:
            short_opt, argdata = copy.deepcopy(argument_data[arg])
        except KeyError:
            print("Skipping undefined pypiper argument: '{}'".format(arg))
            continue
        argdata["required"] = arg in required
        parser.add_argument(short_opt, "--{}".format(arg), **argdata)

    return parser
