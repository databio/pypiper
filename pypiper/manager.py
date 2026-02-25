#!/usr/env python

"""
Pypiper is a Python module for building restartable, robust shell pipelines.

The core component is the PipelineManager class, which wraps shell command
execution with automatic checkpointing, file integrity protection, elapsed time
tracking, and resource monitoring. Domain-specific toolkits (e.g. NGSTk for
genomics) provide optional higher-level convenience functions.
"""

import atexit
import datetime
import errno
import glob
import os
import platform
import re
import shlex  # for splitting commands like a shell does
import signal
import subprocess
import sys
import threading
import time
import warnings
from collections.abc import Callable, Iterable
from hashlib import md5
from importlib.metadata import version
from typing import IO, Any

import logmuse
import pandas as _pd
import psutil
from pipestat import PipestatError, PipestatManager
from yacman import load_yaml

import __main__

from .const import DEFAULT_SAMPLE_NAME, PROFILE_COLNAMES
from .echo_dict import EchoDict

__version__ = version("piper")
__pipestat_version__ = version("pipestat")
from .exceptions import PipelineHalt, SubprocessError
from .flags import *
from .utils import (
    CHECKPOINT_SPECIFICATIONS,
    _check_shell,
    _checkpoint_filepath,
    _clear_flags,
    _default_pipeline_config,
    _flag_name,
    _get_proc_name,
    _is_multi_target,
    _make_lock_name,
    _parse_cmd,
    _pipeline_filepath,
    logger_via_cli,
    result_formatter_markdown,
)

__all__ = ["PipelineManager"]


LOCK_PREFIX = "lock."
LOGFILE_SUFFIX = "_log.md"


class Unbuffered(object):
    def __init__(self, stream: IO) -> None:
        self.stream = stream

    def write(self, data: str) -> None:
        self.stream.write(data)
        self.stream.flush()

    def writelines(self, datas: Iterable[str]) -> None:
        self.stream.writelines(datas)
        self.stream.flush()

    def __getattr__(self, attr: str) -> Any:
        return getattr(self.stream, attr)


class PipelineManager(object):
    """Manage a pipeline: run commands, track files, and report results.

    Example:
        pm = PipelineManager(name="my_pipeline", outfolder="output/")
        pm.run("samtools sort in.bam > out.bam", target="out.bam")
        pm.report_result("read_count", 1000000)
        pm.complete()

    PipelineManager handles command execution with file locking, target-based
    skipping for restartability, intermediate file cleanup, and result reporting
    via pipestat. Each command can declare a target file; if the target exists,
    the command is skipped. Lock files prevent parallel pipelines from colliding
    on the same target. Checkpoints (via timestamp()) enable start/stop control
    for partial reruns.

    Args:
        name: Pipeline name, used for output files, flags, and logging.
        outfolder: Directory for pipeline results.
        version: Pipeline version string.
        args: Parsed argparse.Namespace; pypiper records these and extracts
            relevant options (recover, new_start, etc.).
        multi: Enable multiple pipelines in one script (disables log tee).
        dirty: Never auto-delete intermediate files.
        recover: Overwrite lock files to restart a failed pipeline.
        new_start: Rerun every command even if output exists.
        force_follow: Always run follow functions even if command was skipped.
        cores: Number of processors. Default: 1.
        mem: Memory limit with unit suffix [K|M|G|T]. Default: "1000M".
        config_file: Path to pipeline YAML configuration file.
        output_parent: Parent directory for the output folder.
        overwrite_checkpoints: Ignore checkpoint files (used by Pipeline class).
        logger_kwargs: Keyword arguments for logmuse logger setup.
        pipestat_record_identifier: Record ID for pipestat reporting.
        pipestat_schema: Path to pipestat output schema.
        pipestat_results_file: Path to YAML file backend for results.
        pipestat_config: Path to pipestat configuration file.
        pipestat_pipeline_type: "sample" or "project".
        pipestat_validate_results: Validate results against schema.
            None (default) auto-detects based on schema presence.
        pipestat_additional_properties: Allow results not in schema.
            None (default) uses schema's own setting.
        pipestat_result_formatter: Callable for formatting reported results.
    """

    def __init__(
        self,
        name: str,
        outfolder: str,
        version: str | None = None,
        args: Any = None,
        multi: bool = False,
        dirty: bool = False,
        recover: bool = False,
        new_start: bool = False,
        force_follow: bool = False,
        cores: int = 1,
        mem: str = "1000M",
        config_file: str | None = None,
        output_parent: str | None = None,
        overwrite_checkpoints: bool = False,
        logger_kwargs: dict[str, Any] | None = None,
        pipestat_record_identifier: str | None = None,
        pipestat_schema: str | None = None,
        pipestat_results_file: str | None = None,
        pipestat_config: str | None = None,
        pipestat_pipeline_type: str | None = None,
        pipestat_validate_results: bool | None = None,
        pipestat_additional_properties: bool | None = None,
        pipestat_result_formatter: Callable | None = None,
        **kwargs: Any,
    ) -> None:
        # Params defines the set of options that could be updated via
        # command line args to a pipeline run, that can be forwarded
        # to Pypiper. If any pypiper arguments are passed
        # (via add_pypiper_args()), these will override the constructor
        # defaults for these arguments.

        # Establish default params
        params = {
            "dirty": dirty,
            "recover": recover,
            "new_start": new_start,
            "force_follow": force_follow,
            "config_file": config_file,
            "output_parent": output_parent,
            "cores": cores,
            "mem": mem,
            "testmode": False,
        }

        # Transform the command-line namespace into a Mapping.
        args_dict = vars(args) if args else dict()

        # Parse and store stage specifications that can determine pipeline
        # start and/or stop point.
        # First, add such specifications to the command-line namespace,
        # favoring the command-line spec if both are present.
        for cp_spec in set(CHECKPOINT_SPECIFICATIONS) & set(kwargs.keys()):
            args_dict.setdefault(cp_spec, kwargs[cp_spec])
        # Then, ensure that we set each such specification on this manager
        # so that we're guaranteed safe attribute access. If it's present,
        # remove the specification from the namespace that will be used to
        # update this manager's parameters Mapping.
        for optname in CHECKPOINT_SPECIFICATIONS:
            checkpoint = args_dict.pop(optname, None)
            setattr(self, optname, checkpoint)
        if self.stop_before and self.stop_after:
            raise TypeError(
                "Cannot specify both stop_before and stop_after. "
                "stop_before='{before}' (exclusive: stage is NOT run) and "
                "stop_after='{after}' (inclusive: stage IS run) are mutually exclusive. "
                "Use only one.".format(before=self.stop_before, after=self.stop_after)
            )

        # Update this manager's parameters with non-checkpoint-related
        # command-line parameterization.
        params.update(args_dict)

        # If no starting point was specified, assume that the pipeline's
        # execution is to begin right away and set the internal flag so that
        # run() is let loose to execute instructions given.
        self._active = not self.start_point
        self._stopped = False

        # Pipeline-level variables to track global state and pipeline stats
        # Pipeline settings
        self.name = name
        self.tee = None
        self.overwrite_locks = params["recover"]
        self.new_start = params["new_start"]
        self.force_follow = params["force_follow"]
        self.dirty = params["dirty"]
        self.cores = params["cores"]
        self.output_parent = params["output_parent"]
        self.testmode = params["testmode"]

        # Establish the log file to check safety with logging keyword arguments.
        # Establish the output folder since it's required for the log file.
        self.outfolder = os.path.join(outfolder, "")  # trailing slash
        self.pipeline_log_file = _pipeline_filepath(self, suffix=LOGFILE_SUFFIX)

        # Set up logger
        logger_kwargs = logger_kwargs or {}
        if logger_kwargs.get("logfile") == self.pipeline_log_file:
            raise ValueError(
                f"The logfile given for the pipeline manager's logger matches that which will be used by the manager itself: {self.pipeline_log_file}"
            )
        default_logname = ".".join([__name__, self.__class__.__name__, self.name])
        self._logger = None
        if args:
            logger_builder_method = "logger_via_cli"
            try:
                self._logger = logger_via_cli(args, **logger_kwargs)
            except logmuse.est.AbsentOptionException as e:
                # Defer logger construction to init_logger.
                self.debug(f"logger_via_cli failed: {e}")
        if self._logger is None:
            logger_builder_method = "init_logger"
            # covers cases of bool(args) being False, or failure of logger_via_cli.
            # strict is only for logger_via_cli.
            logger_kwargs = {k: v for k, v in logger_kwargs.items() if k != "strict"}
            try:
                name = logger_kwargs.pop("name")
            except KeyError:
                name = default_logname
            self._logger = logmuse.init_logger(name, **logger_kwargs)
        self.debug(f"Logger set with {logger_builder_method}")

        # Keep track of an ID for the number of processes attempted
        self.proc_count = 0

        # We use this memory to pass a memory limit to processes like java that
        # can take a memory limit, so they don't get killed by a SLURM (or other
        # cluster manager) overage. However, with java, the -Xmx argument can only
        # limit the *heap* space, not total memory use; so occasionally SLURM will
        # still kill these processes because total memory goes over the limit.
        # As a kind of hack, we'll set the java processes heap limit to 95% of the
        # total memory limit provided.
        # This will give a little breathing room for non-heap java memory use.

        if not params["mem"].endswith(("K", "M", "G", "T")):
            self.mem = params["mem"] + "M"
        else:
            # Assume the memory is in megabytes.
            self.mem = params["mem"]

        self.javamem = str(int(int(self.mem[:-1]) * 0.95)) + self.mem[-1:]

        self.container = None
        self.clean_initialized = False

        # Do some cores math for split processes
        # If a pipeline wants to run a process using half the cores, or 1/4 of the cores,
        # this can lead to complications if the number of cores is not evenly divisible.
        # Here we add a few variables so that pipelines can easily divide the cores evenly.
        # 50/50 split
        self.cores1of2a = int(self.cores) / 2 + int(self.cores) % 2
        self.cores1of2 = int(self.cores) / 2

        # 75/25 split
        self.cores1of4 = int(self.cores) / 4
        self.cores3of4 = int(self.cores) - int(self.cores1of4)

        self.cores1of8 = int(self.cores) / 8
        self.cores7of8 = int(self.cores) - int(self.cores1of8)

        self.pl_version = version
        # Set relative output_parent directory to absolute
        # not necessary after all. . .
        # if self.output_parent and not os.path.isabs(self.output_parent):
        #   self.output_parent = os.path.join(os.getcwd(), self.output_parent)

        # File paths:
        self.make_sure_path_exists(self.outfolder)
        self.pipeline_profile_file = _pipeline_filepath(self, suffix="_profile.tsv")

        # Stats and figures are general and so lack the pipeline name.
        self.pipeline_stats_file = _pipeline_filepath(self, filename="stats.yaml")

        # Record commands used and provide manual cleanup script.
        self.pipeline_commands_file = _pipeline_filepath(self, suffix="_commands.sh")
        self.cleanup_file = _pipeline_filepath(self, suffix="_cleanup.sh")

        # Pipeline status variables
        self.peak_memory = 0  # memory high water mark
        self.starttime = time.time()
        self.last_timestamp = self.starttime  # time of the last call to timestamp()

        self.locks = []
        self.running_procs = {}
        self.completed_procs = {}

        self.wait = True  # turn off for debugging

        # Initialize status and flags
        self.status = "initializing"
        # as part of the beginning of the pipeline, clear any flags set by
        # previous runs of this pipeline
        _clear_flags(self)

        # In-memory holder for report_result
        self.stats_dict = {}

        # Result formatter to pass to pipestat
        self.pipestat_result_formatter = pipestat_result_formatter or result_formatter_markdown

        # Checkpoint-related parameters
        self.overwrite_checkpoints = overwrite_checkpoints or self.new_start
        self.halt_on_next = False
        self.prev_checkpoint = None
        self.curr_checkpoint = None

        # Pypiper can keep track of intermediate files to clean up at the end
        self.cleanup_list = []
        self.cleanup_list_conditional = []

        # Register handler functions to deal with interrupt and termination signals;
        # If received, we would then clean up properly (set pipeline status to FAIL, etc).
        # signal.signal() only works in the main thread; skip if called from a
        # worker thread (e.g. AI agent tool calls, web servers, thread pools).
        if threading.current_thread() is threading.main_thread():
            signal.signal(signal.SIGINT, self._signal_int_handler)
            signal.signal(signal.SIGTERM, self._signal_term_handler)

        # pipestat setup
        self.pipestat_record_identifier = pipestat_record_identifier or DEFAULT_SAMPLE_NAME
        self.pipestat_pipeline_type = pipestat_pipeline_type or "sample"

        # don't force default pipestat_results_file value unless
        # pipestat config not provided
        if pipestat_config is None and pipestat_results_file is None:
            self.pipestat_results_file = self.pipeline_stats_file
        elif pipestat_results_file:
            self.pipestat_results_file = pipestat_results_file
            self.pipeline_stats_file = self.pipestat_results_file

        def _get_arg(args_dict, arg_name):
            """safely get argument from arg dict -- return None if doesn't exist"""
            return None if arg_name not in args_dict else args_dict[arg_name]

        # Resolve the schema path from explicit arg or CLI
        resolved_schema = pipestat_schema or _get_arg(args_dict, "pipestat_schema")

        # Resolve pipestat_validate_results from constructor arg or CLI
        if pipestat_validate_results is None:
            cli_val = _get_arg(args_dict, "pipestat_validate_results")
            if cli_val is not None:
                pipestat_validate_results = cli_val.lower() == "true"

        # Resolve pipestat_additional_properties from constructor arg or CLI
        if pipestat_additional_properties is None:
            cli_val = _get_arg(args_dict, "pipestat_additional_properties")
            if cli_val is not None:
                pipestat_additional_properties = cli_val.lower() == "true"

        # Resolve validate_results: if explicitly provided, use that value.
        # Otherwise, pass None to let PipestatManager auto-detect from schema presence.
        resolved_validate = pipestat_validate_results

        pipestat_config_resolved = pipestat_config or _get_arg(args_dict, "pipestat_config")
        if pipestat_config_resolved:
            self._pipestat_manager = PipestatManager.from_config(
                config=pipestat_config_resolved,
                pipeline_type=self.pipestat_pipeline_type,
                multi_pipelines=multi,
            )
        else:
            self._pipestat_manager = PipestatManager.from_file_backend(
                results_file_path=self.pipestat_results_file
                or _get_arg(args_dict, "pipestat_results_file")
                or self.pipeline_stats_file,
                schema_path=resolved_schema,
                record_identifier=self.pipestat_record_identifier
                or _get_arg(args_dict, "pipestat_sample_name")
                or DEFAULT_SAMPLE_NAME,
                pipeline_name=self.name,
                pipeline_type=self.pipestat_pipeline_type,
                multi_pipelines=multi,
                validate_results=resolved_validate,
                additional_properties=pipestat_additional_properties,
            )

        # Set result formatter as property (removed from __init__)
        if pipestat_result_formatter:
            self._pipestat_manager.result_formatter = pipestat_result_formatter

        self.start_pipeline(args, multi)

        # Handle config file if it exists

        # Read YAML config file
        # TODO: This section should become a function, so toolkits can use it
        # to locate a config file.
        config_to_load = None  # start with nothing

        if config_file:
            config_to_load = config_file
        else:
            cmdl_config_file = getattr(args, "config_file", None)
            if cmdl_config_file:
                if os.path.isabs(cmdl_config_file):
                    # Absolute custom config file specified
                    if os.path.isfile(cmdl_config_file):
                        config_to_load = cmdl_config_file
                    else:
                        self.debug("Can't find custom config file: " + cmdl_config_file)
                        pass
                else:
                    # Relative custom config file specified
                    # Set path to be relative to pipeline script
                    pipedir = os.path.dirname(sys.argv[0])
                    abs_config = os.path.join(pipedir, cmdl_config_file)
                    if os.path.isfile(abs_config):
                        config_to_load = abs_config
                    else:
                        self.debug("File: {}".format(__file__))
                        self.debug("Can't find custom config file: " + abs_config)
                        pass
                if config_to_load is not None:
                    pass
                    self.debug("\nUsing custom config file: {}".format(config_to_load))
            else:
                # No custom config file specified. Check for default
                default_config = _default_pipeline_config(sys.argv[0])
                if os.path.isfile(default_config):
                    config_to_load = default_config
                    self.debug("Using default pipeline config file: {}".format(config_to_load))

        # Finally load the config we found.
        if config_to_load is not None:
            self.debug("\nLoading config file: {}\n".format(config_to_load))
            self.config = EchoDict(load_yaml(config_to_load))
            # Ensure standard sections exist for nested attribute access
            # (EchoDict returns strings for missing keys, which breaks setting)
            self.config.setdefault("tools", {})
            self.config.setdefault("parameters", {})
            self.config.setdefault("resources", {})
        else:
            self.debug("No config file")
            self.config = None

    @property
    def pipestat(self) -> PipestatManager:
        """Access the PipestatManager for reporting results and managing status.

        Example:
            pm.pipestat.report(values={"reads": 1000}, record_identifier="sample1")
            status = pm.pipestat.get_status("sample1")

        Returns:
            Configured PipestatManager instance.

        Raises:
            PipestatError: If pipestat was not initialized (no schema provided).
        """
        try:
            return getattr(self, "_pipestat_manager")
        except AttributeError:
            raise PipestatError(
                "PipestatManager has not been configured for this pipeline run. "
                "To enable pipestat result reporting, provide pipestat_schema='/path/to/output_schema.yaml' "
                "when creating the PipelineManager, or pass --pipestat-schema on the command line."
            )

    @property
    def _completed(self) -> bool:
        """Whether the managed pipeline is in a completed state."""
        return self.pipestat.get_status(self._pipestat_manager.record_identifier) == COMPLETE_FLAG

    @property
    def _failed(self) -> bool:
        """Whether the managed pipeline is in a failed state."""
        return self.pipestat.get_status(self._pipestat_manager.record_identifier) == FAIL_FLAG

    @property
    def halted(self) -> bool:
        """Whether the managed pipeline is in a paused/halted state."""
        return self.pipestat.get_status(self._pipestat_manager.record_identifier) == PAUSE_FLAG

    @property
    def _has_exit_status(self) -> bool:
        """Whether the managed pipeline has been safely stopped."""
        if self._stopped:
            return True
        return self._completed or self.halted or self._failed

    def _ignore_interrupts(self) -> None:
        """Ignore interrupt and termination signals for subprocess preexec_fn."""
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGTERM, signal.SIG_IGN)

    def start_pipeline(self, args: Any = None, multi: bool = False) -> None:
        """Initialize pipeline logging, diagnostics, and status tracking.

        Called automatically by __init__; rarely called directly.

        Sets up log file tee (mirroring stdout to a log file), prints
        version/environment info, records git state, creates output folder,
        and sets pipeline status to 'running'. In interactive/multi mode,
        the log tee is disabled to avoid interfering with notebook or
        multi-pipeline output.
        """
        # Perhaps this could all just be put into __init__, but I just kind of like the idea of a start function
        # self.make_sure_path_exists(self.outfolder)

        # By default, Pypiper will mirror every operation so it is displayed both
        # on sys.stdout **and** to a log file. Unfortunately, interactive python sessions
        # ruin this by interfering with stdout. So, for interactive mode, we do not enable
        # the tee subprocess, sending all output to screen only.
        # Starting multiple PipelineManagers in the same script has the same problem, and
        # must therefore be run in interactive_mode.

        interactive_mode = multi or not hasattr(__main__, "__file__")
        if interactive_mode:
            self.warning(
                "Warning: You're running an interactive python session. "
                "This works, but pypiper cannot tee the output, so results "
                "are only logged to screen."
            )
        else:
            sys.stdout = Unbuffered(sys.stdout)
            # sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)  # Unbuffer output

            # The tee subprocess must be instructed to ignore TERM and INT signals;
            # Instead, I will clean up this process in the signal handler functions.
            # This is required because otherwise, if pypiper receives a TERM or INT,
            # the tee will be automatically terminated by python before I have a chance to
            # print some final output (for example, about when the process stopped),
            # and so those things don't end up in the log files because the tee
            # subprocess is dead. Instead, I will handle the killing of the tee process
            # manually (in the exit handler).

            # a for append to file

            tee = subprocess.Popen(
                ["tee", "-a", self.pipeline_log_file],
                stdin=subprocess.PIPE,
                preexec_fn=self._ignore_interrupts,
            )

            # If the pipeline is terminated with SIGTERM/SIGINT,
            # make sure we kill this spawned tee subprocess as well.
            # atexit.register(self._kill_child_process, tee.pid, proc_name="tee")
            os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
            os.dup2(tee.stdin.fileno(), sys.stderr.fileno())

            self.tee = tee

        # For some reason, this exit handler function MUST be registered after
        # the one that kills the tee process.
        atexit.register(self._exit_handler)

        # A future possibility to avoid this tee, is to use a Tee class; this works for anything printed here
        # by pypiper, but can't tee the subprocess output. For this, it would require using threading to
        # simultaneously capture and display subprocess output. I shelve this for now and stick with the tee option.
        # sys.stdout = Tee(self.pipeline_log_file)

        # Record the git version of the pipeline and pypiper used. This gets (if it is in a git repo):
        # dir: the directory where the code is stored
        # hash: the commit id of the last commit in this repo
        # date: the date of the last commit in this repo
        # diff: a summary of any differences in the current (run) version vs. the committed version

        # Wrapped in try blocks so that the code will not fail if the pipeline or pypiper are not git repositories
        gitvars = {}
        try:
            # pypiper dir
            ppd = os.path.dirname(os.path.realpath(__file__))
            gitvars["pypiper_dir"] = ppd
            gitvars["pypiper_hash"] = (
                subprocess.check_output(
                    "cd " + ppd + "; git rev-parse --verify HEAD 2>/dev/null",
                    shell=True,
                )
                .decode()
                .strip()
            )
            gitvars["pypiper_date"] = (
                subprocess.check_output(
                    "cd " + ppd + "; git show -s --format=%ai HEAD 2>/dev/null",
                    shell=True,
                )
                .decode()
                .strip()
            )
            gitvars["pypiper_diff"] = (
                subprocess.check_output(
                    "cd " + ppd + "; git diff --shortstat HEAD 2>/dev/null", shell=True
                )
                .decode()
                .strip()
            )
            gitvars["pypiper_branch"] = (
                subprocess.check_output(
                    "cd " + ppd + "; git branch | grep '*' 2>/dev/null", shell=True
                )
                .decode()
                .strip()
            )
        except Exception:
            pass
        try:
            # pipeline dir
            pld = os.path.dirname(os.path.realpath(sys.argv[0]))
            gitvars["pipe_dir"] = pld
            gitvars["pipe_hash"] = (
                subprocess.check_output(
                    "cd " + pld + "; git rev-parse --verify HEAD 2>/dev/null",
                    shell=True,
                )
                .decode()
                .strip()
            )
            gitvars["pipe_date"] = (
                subprocess.check_output(
                    "cd " + pld + "; git show -s --format=%ai HEAD 2>/dev/null",
                    shell=True,
                )
                .decode()
                .strip()
            )
            gitvars["pipe_diff"] = (
                subprocess.check_output(
                    "cd " + pld + "; git diff --shortstat HEAD 2>/dev/null", shell=True
                )
                .decode()
                .strip()
            )
            gitvars["pipe_branch"] = (
                subprocess.check_output(
                    "cd " + pld + "; git branch | grep '*' 2>/dev/null", shell=True
                )
                .decode()
                .strip()
            )
        except Exception:
            pass

        # Print out a header section in the pipeline log:
        # Wrap things in backticks to prevent markdown from interpreting underscores as emphasis.
        # print("----------------------------------------")
        def logfmt(key, value=None, padding=16):
            padded_key = key.rjust(padding)
            formatted_val = f"`{value}`" if value else ""
            return f"* {padded_key}: {formatted_val}"

        self.info("### Pipeline run code and environment:\n")
        self.info(logfmt("Command", str(" ".join(sys.argv))))
        self.info(logfmt("Compute host", platform.node()))
        self.info(logfmt("Working dir", os.getcwd()))
        self.info(logfmt("Outfolder", self.outfolder))
        self.info(logfmt("Log file", self.pipeline_log_file))
        self.timestamp(logfmt("Start time"))

        self.info("\n### Version log:\n")
        self.info(logfmt("Python version", platform.python_version()))
        try:
            self.info(logfmt("Pypiper dir", gitvars["pypiper_dir"].strip()))
            self.info(logfmt("Pypiper version", __version__))
            self.info(logfmt("Pypiper hash", gitvars["pypiper_hash"]))
            self.info(logfmt("Pypiper branch", gitvars["pypiper_branch"]))
            self.info(logfmt("Pypiper date", gitvars["pypiper_date"]))
            if gitvars["pypiper_diff"]:
                self.info(logfmt("Pypiper diff", gitvars["pypiper_diff"]))
        except KeyError:
            # It is ok if keys aren't set, it means pypiper isn't in a  git repo.
            pass

        self.info(logfmt("Pipestat version", __pipestat_version__))

        try:
            self.info(logfmt("Pipeline dir", gitvars["pipe_dir"].strip()))
            self.info(logfmt("Pipeline version", self.pl_version))
            self.info(logfmt("Pipeline hash", gitvars["pipe_hash"]).strip())
            self.info(logfmt("Pipeline branch", gitvars["pipe_branch"]).strip())
            self.info(logfmt("Pipeline date", gitvars["pipe_date"]).strip())
            if gitvars["pipe_diff"] != "":
                self.info(logfmt("Pipeline diff", gitvars["pipe_diff"]).strip())
        except KeyError:
            # It is ok if keys aren't set, it means the pipeline isn't a git repo.
            pass

        # self.info all arguments (if any)
        self.info("\n### Arguments passed to pipeline:\n")
        for arg, val in sorted((vars(args) if args else dict()).items()):
            argtext = "`{}`".format(arg)
            valtext = "`{}`".format(val)
            self.info("* {}:  {}".format(argtext.rjust(20), valtext))

        self.info("\n### Initialized Pipestat Object:\n")
        results = self._pipestat_manager.__str__().split("\n")
        for i in results:
            self.info("* " + i)
        self.info("* Sample name: " + self.pipestat_record_identifier + "\n")
        self.info("\n----------------------------------------\n")
        self.status = "running"
        self.pipestat.set_status(
            record_identifier=self._pipestat_manager.record_identifier,
            status_identifier="running",
        )

        # Record the start in PIPE_profile and PIPE_commands output files so we
        # can trace which run they belong to

        with open(self.pipeline_commands_file, "a") as myfile:
            myfile.write(
                "# Pipeline started at "
                + time.strftime("%m-%d %H:%M:%S", time.localtime(self.starttime))
                + "\n\n"
            )

        with open(self.pipeline_profile_file, "a") as myfile:
            myfile.write(
                "# Pipeline started at "
                + time.strftime("%m-%d %H:%M:%S", time.localtime(self.starttime))
                + "\n\n"
                + "# "
                + "\t".join(PROFILE_COLNAMES)
                + "\n"
            )

    def _set_status_flag(self, status: str) -> None:
        """
        Configure state and files on disk to match current processing status.

        Args:
            status (str): Name of new status designation for pipeline.
        """

        # Remove previous status flag file.
        flag_file_path = self._flag_file_path()
        try:
            os.remove(flag_file_path)
        except Exception:
            # Print message only if the failure to remove the status flag
            # is unexpected; there's no flag for initialization, so we
            # can't remove the file.
            if self.status != "initializing":
                self.debug("Could not remove flag file: '{}'".format(flag_file_path))
            pass

        # Set new status.
        prev_status = self.status
        self.status = status
        self.pipestat.set_status(
            record_identifier=self._pipestat_manager.record_identifier,
            status_identifier=status,
        )
        self.debug("\nChanged status from {} to {}.".format(prev_status, self.status))

    def _flag_file_path(self, status: str | None = None) -> str:
        """
        Create path to flag file based on indicated or current status.

        Internal variables used are the pipeline name and the designated
        pipeline output folder path.

        Args:
            status (str): flag file type to create, default to current status

        Returns:
            str: path to flag file of indicated or current status.
        """

        flag_file_name = "{}_{}_{}".format(
            self._pipestat_manager.pipeline_name,
            self.pipestat_record_identifier,
            _flag_name(status or self.status),
        )
        return _pipeline_filepath(self, filename=flag_file_name)

    ###################################
    # Process calling functions
    ###################################
    def run(
        self,
        cmd: str | list[str],
        target: str | list[str] | None = None,
        lock_name: str | list[str] | None = None,
        shell: bool | None = None,
        nofail: bool = False,
        clean: bool = False,
        follow: Callable | None = None,
        container: str | None = None,
        default_return_code: int | None = 0,
    ) -> int | None:
        """Run a shell command with file-locking and target-based skipping.

        Example:
            pm.run("sort input.txt > output.txt", target="output.txt")
            pm.run(["cmd1", "cmd2"], target="final.out")
            pm.run("samtools index in.bam", target="in.bam.bai", clean=True)

        If target exists, the command is skipped (restartability). If another
        process holds the lock, waits for it to finish. Creates lock files
        during execution to prevent parallel conflicts.

        Follow functions run only when the command actually executes (target
        didn't exist), unless force_follow=True was set on the manager.

        Args:
            cmd: Shell command string, or list of commands to run in sequence.
                A list of commands runs each sequentially, returning the max
                return code.
            target: Output file(s). If all exist, command is skipped.
                If None, lock_name is required. Multiple targets can be
                provided as a list; all must exist to skip.
            lock_name: Explicit lock file name. Defaults to target-based name.
                Required if target is None (for commands with no output file).
            shell: Force shell mode. None (default) auto-detects based on
                the presence of pipes (|) or redirects (>).
            nofail: If True, pipeline continues past nonzero return codes.
                The failure is logged but does not halt the pipeline.
            clean: If True, adds target to auto-cleanup list (deleted on
                pipeline success, kept with --dirty).
            follow: Callable to run after command execution. Skipped if the
                command was skipped (target existed), unless force_follow=True.
            container: Docker container name for execution.
            default_return_code: Return code when command is skipped (target
                existed). Default: 0.

        Returns:
            Return code (int). For command lists, the maximum return code.
        """

        def _max_ret_code(codes_list):
            """
            Return the maximum of a list of return codes.

            Args:
                code (list[int]): List of return codes to compare.

            Returns:
                int: Maximum of list.
            """
            # filter out codes that are None
            codes_list = [code for code in codes_list if code is not None]
            # get the max of the remaining codes
            if codes_list:
                return max(codes_list)
            # if no codes are left, return None
            return

        # validate default return code
        if default_return_code is not None and not isinstance(default_return_code, int):
            raise TypeError("default_return_code must be an int or None")

        # If the pipeline's not been started, skip ahead.
        if not self._active:
            cmds = [cmd] if isinstance(cmd, str) else cmd
            cmds_text = [c if isinstance(c, str) else " ".join(c) for c in cmds]
            self.info(
                "Pipeline is inactive; skipping {} command(s):\n{}".format(
                    len(cmds), "\n".join(cmds_text)
                )
            )
            return default_return_code

        # Short-circuit if the checkpoint file exists and the manager's not
        # been configured to overwrite such files.
        if self.curr_checkpoint is not None:
            check_fpath = _checkpoint_filepath(self.curr_checkpoint, self)
            if os.path.isfile(check_fpath) and not self.overwrite_checkpoints:
                self.info(
                    "Checkpoint file exists for '{stage}' at '{path}'. "
                    "Skipping command: `{cmd}`. "
                    "To re-run this stage, delete the checkpoint file or use new_start=True (CLI: -N).".format(
                        stage=self.curr_checkpoint, path=check_fpath, cmd=cmd
                    )
                )
                return default_return_code

        # TODO: consider making the logic such that locking isn't implied, or
        # TODO (cont.): that we can make it otherwise such that it's not
        # TODO (cont.): strictly necessary to provide target or lock_name.
        # The default lock name is based on the target name.
        # Therefore, a targetless command that you want
        # to lock must specify a lock_name manually.
        if target is None and lock_name is None:
            self.fail_pipeline(
                Exception(
                    "PipelineManager.run() requires either a 'target' (output file path) or a 'lock_name'. "
                    "Provide target='/path/to/output_file' to enable output checking and file locking, "
                    "or provide lock_name='my_step' for targetless commands that still need locking."
                )
            )

        # Downstream code requires target to be a list, so convert if only
        # a single item was given
        if not _is_multi_target(target) and target is not None:
            target = [target]

        # Downstream code requires a list of locks; convert
        if isinstance(lock_name, str):
            lock_name = [lock_name]

        # Default lock_name (if not provided) is based on the target file name,
        # but placed in the parent pipeline outfolder
        self.debug(
            "Lock_name {}; target '{}', outfolder '{}'".format(lock_name, target, self.outfolder)
        )
        lock_name = lock_name or _make_lock_name(target, self.outfolder)
        lock_files = [self._make_lock_path(ln) for ln in lock_name]

        process_return_code = default_return_code
        local_maxmem = 0

        # Decide how to do follow-up.
        if not follow:

            def call_follow():
                return None
        elif not hasattr(follow, "__call__"):
            # Warn about non-callable argument to follow-up function.
            self.warning(
                "Follow-up function is not callable and won't be used: {}".format(type(follow))
            )

            def call_follow():
                return None
        else:
            # Wrap the follow-up function so that the log shows what's going on.
            # additionally, the in_follow attribute is set to enable proper command count handling
            def call_follow():
                self.debug("Follow:")
                self.in_follow = True
                follow()
                self.in_follow = False

        # The while=True loop here is unlikely to be triggered, and is just a
        # wrapper to prevent race conditions; the lock_file must be created by
        # the current loop. If not, we loop again and then re-do the tests.
        # The recover and newstart options inform the pipeline to run a command
        # in a scenario where it normally would not. We use these "local" flags
        # to allow us to report on the state of the pipeline in the first round
        # as normal, but then proceed on the next iteration through the outer
        # loop. The proceed_through_locks is a flag that is set if any lockfile
        # is found that needs to be recovered or overwritten. It instructs us to
        # ignore lock files on the next iteration.
        local_recover = False
        local_newstart = False
        proceed_through_locks = False

        while True:
            ##### Tests block
            # Base case: All targets exists and not set to overwrite targets break loop, don't run process.
            # os.path.exists returns True for either a file or directory; .isfile is file-only
            if (
                target is not None
                and all([os.path.exists(t) for t in target])
                and not any([os.path.isfile(lf) for lf in lock_files])
                and not local_newstart
            ):
                for tgt in target:
                    if os.path.exists(tgt):
                        self.info(
                            "Target exists: `{tgt}`. Skipping this step. "
                            "To force re-computation, use new_start=True (CLI: -N).".format(
                                tgt=tgt
                            )
                        )
                if self.new_start:
                    self.info("New start mode; run anyway.  ")
                    # Set the local_newstart flag so the command will run anyway.
                    # Doing this in here instead of outside the loop allows us
                    # to still report the target existence.
                    local_newstart = True
                    continue
                # Normally we don't run the follow, but if you want to force. . .
                if self.force_follow:
                    call_follow()
                # Increment process count
                increment_info_pattern = (
                    "Skipped command: `{}`\nCommand ID incremented by: `{}`. Current ID: `{}`\n"
                )
                if isinstance(cmd, list):
                    for c in cmd:
                        count = len(_parse_cmd(c, shell))
                        self.proc_count += count
                        self.debug(increment_info_pattern.format(str(c), count, self.proc_count))
                else:
                    count = len(_parse_cmd(cmd, shell))
                    self.proc_count += count
                    self.debug(increment_info_pattern.format(str(cmd), count, self.proc_count))
                break  # Do not run command

            # Scenario 1: Lock file exists, but we're supposed to overwrite target; Run process.
            if not proceed_through_locks:
                for lock_file in lock_files:
                    recover_file = self._recoverfile_from_lockfile(lock_file)
                    if os.path.isfile(lock_file):
                        self.info(
                            "Found lock file: {lock}. "
                            "This means another pipeline may be running on this target, "
                            "or a previous run crashed without cleaning up. "
                            "To override the lock and re-run the command (overwriting any partial output), "
                            "restart with recover=True (CLI: -R).".format(lock=lock_file)
                        )
                        if self.overwrite_locks:
                            self.info("Overwriting target...")
                            proceed_through_locks = True
                        elif os.path.isfile(recover_file):
                            self.info(
                                "Found dynamic recovery file ({}); overwriting target...".format(
                                    recover_file
                                )
                            )
                            # remove the lock file which will then be promptly re-created for the current run.
                            local_recover = True
                            proceed_through_locks = True
                            # the recovery flag is now spent; remove so we don't accidentally re-recover a failed job
                            os.remove(recover_file)
                        else:  # don't overwrite locks
                            self._wait_for_lock(lock_file)
                            # when it's done loop through again to try one more
                            # time (to see if the target exists now)
                            continue

            # If you get to this point, the target doesn't exist, and the lock_file doesn't exist
            # (or we should overwrite). create the lock (if you can)
            # Initialize lock in master lock list
            for lock_file in lock_files:
                self.locks.append(lock_file)
                if self.overwrite_locks or local_recover:
                    self._create_file(lock_file)
                else:
                    try:
                        self._create_file_racefree(lock_file)  # Create lock
                    except OSError as e:
                        if e.errno == errno.EEXIST:  # File already exists
                            self.info(
                                "Lock file appeared between existence check and creation (race condition): {lock}. "
                                "Re-checking. This is normal when multiple pipelines target the same file.".format(
                                    lock=lock_file
                                )
                            )

                            # Since a lock file was created by a different source,
                            # we need to reset this flag to re-check the locks.
                            proceed_through_locks = False
                            continue  # Go back to start

            ##### End tests block
            # If you make it past these tests, we should proceed to run the process.

            if target is not None:
                self.info(
                    "Target to produce: {}  ".format(",".join(["`" + x + "`" for x in target]))
                )
            else:
                self.info("Targetless command, running...  ")

            if isinstance(cmd, list):  # Handle command lists
                for cmd_i in cmd:
                    list_ret, maxmem = self.callprint(cmd_i, shell, lock_file, nofail, container)
                    maxmem = max(maxmem) if isinstance(maxmem, Iterable) else maxmem
                    local_maxmem = max(local_maxmem, maxmem)
                    list_ret = (
                        _max_ret_code(list_ret) if isinstance(list_ret, Iterable) else list_ret
                    )
                    process_return_code = _max_ret_code([process_return_code, list_ret])

            else:  # Single command (most common)
                process_return_code, local_maxmem = self.callprint(
                    cmd, shell, lock_file, nofail, container
                )  # Run command
                if isinstance(process_return_code, list):
                    process_return_code = _max_ret_code(process_return_code)

            # For temporary files, you can specify a clean option to automatically
            # add them to the clean list, saving you a manual call to clean_add
            if target is not None and clean:
                for tgt in target:
                    self.clean_add(tgt)

            call_follow()
            for lock_file in lock_files:
                os.remove(lock_file)  # Remove lock file
                self.locks.remove(lock_file)

            # If you make it to the end of the while loop, you're done
            break

        return process_return_code

    def checkprint(self, cmd: str, shell: bool | None = None, nofail: bool = False) -> str:
        """Run a command and return its stdout as a string.

        Example:
            genome_size = pm.checkprint("wc -l genome.fa")
            version = pm.checkprint("samtools --version")

        Like callprint, but captures and returns stdout (uses
        subprocess.check_output instead of subprocess.call). Use this when
        you need a command's output as a Python variable.

        Args:
            cmd: Shell command string.
            shell: Force shell mode. None auto-detects based on pipes/redirects.
            nofail: If True, pipeline continues past failure instead of halting.

        Returns:
            Stripped stdout string from the command. Empty string in test mode.
        """

        self._report_command(cmd)
        if self.testmode:
            return ""

        likely_shell = _check_shell(cmd, shell)

        if shell is None:
            shell = likely_shell

        if not shell:
            if likely_shell:
                self.debug(
                    "Should this command run in a shell instead of directly in a subprocess?"
                )
            cmd = shlex.split(cmd)

        try:
            return subprocess.check_output(cmd, shell=shell).decode().strip()
        except Exception as e:
            self._triage_error(e, nofail)

    def _attend_process(self, proc: Any, sleeptime: float) -> bool:
        """
        Waits on a process for a given time to see if it finishes, returns True
        if it's still running after the given time or False as soon as it
        returns.

        Args:
            proc (psutil.Popen): Process object opened by psutil.Popen()
            sleeptime (float): Time to wait

        Returns:
            bool: True if process is still running; otherwise false
        """
        # print("attend:{}".format(proc.pid))
        try:
            proc.wait(timeout=sleeptime)
        except psutil.TimeoutExpired:
            return True
        return False

    def callprint(
        self,
        cmd: str,
        shell: bool | None = None,
        lock_file: str | None = None,
        nofail: bool = False,
        container: str | None = None,
    ) -> tuple[list[int | None], list[float]]:
        """Execute a command, print it, and track memory usage.

        Example:
            retcode, maxmem = pm.callprint("samtools sort in.bam")

        Piped commands (|) are split into subprocesses so each can be memory-
        profiled independently. Shell mode (used for redirects > or wildcards *)
        runs the whole command in a shell, which prevents per-process memory
        tracking. Prefer shell=False (default) when possible for better
        profiling.

        This is the low-level execution method; most users should use run()
        instead, which adds target-based skipping and file locking on top.

        Args:
            cmd: Shell command string.
            shell: Force shell mode. None auto-detects from pipes/redirects.
            lock_file: Lock file name for this execution.
            nofail: If True, pipeline continues past nonzero return codes.
            container: Docker container name for execution.

        Returns:
            Tuple of (return_codes_per_process, peak_memory_GB_per_process).
        """
        # The Popen shell argument works like this:
        # if shell=False, then we format the command (with split()) to be a list of command and its arguments.
        # Split the command to use shell=False;
        # leave it together to use shell=True;

        def get_mem_child_sum(proc):
            try:
                # get children processes
                children = proc.children(recursive=True)
                # get RSS memory of each child proc and sum all
                mem_sum = proc.memory_info().rss
                if children:
                    mem_sum += sum([x.memory_info().rss for x in children])
                # return in gigs
                return mem_sum / 1e9
            except (psutil.NoSuchProcess, psutil.ZombieProcess) as e:
                self.warning(e)
                self.warning("Warning: couldn't add memory use for process: {}".format(proc.pid))
                return 0

        def display_memory(memval):
            return None if memval < 0 else "{}GB".format(round(memval, 3))

        def make_hash(o):
            """
            Convert the object to string and hash it, return None in case of failure

            Args:
                o: object of any type, in our case it is a dict

            Returns:
                str: hashed string representation of the dict
            """
            try:
                hsh = md5(str(o).encode("utf-8")).hexdigest()[:10]
            except Exception as e:
                self.debug(
                    "Could not create hash for '{}', caught exception: {}".format(
                        str(o), e.__class__.__name__
                    )
                )
                hsh = None
            return hsh

        if container:
            cmd = "docker exec " + container + " " + cmd

        if self.testmode:
            self._report_command(cmd)
            return 0, 0

        self.debug("Command: {}".format(cmd))
        param_list = _parse_cmd(cmd, shell)
        # cast all commands to str and concatenate for hashing
        conc_cmd = "".join([str(x["args"]) for x in param_list])
        self.debug("Hashed command '{}': {}".format(conc_cmd, make_hash(conc_cmd)))
        processes = []
        running_processes = []
        completed_processes = []
        start_time = time.time()
        for i in range(len(param_list)):
            running_processes.append(i)
            if i == 0:
                processes.append(psutil.Popen(preexec_fn=os.setsid, **param_list[i]))
            else:
                param_list[i]["stdin"] = processes[i - 1].stdout
                processes.append(psutil.Popen(preexec_fn=os.setsid, **param_list[i]))
            self.running_procs[processes[-1].pid] = {
                "proc_name": _get_proc_name(param_list[i]["args"]),
                "start_time": start_time,
                "container": container,
                "p": processes[-1],
                "args_hash": make_hash(conc_cmd),
                "local_proc_id": self.process_counter(),
            }

        self._report_command(cmd, [x.pid for x in processes])
        # Capture the subprocess output in <pre> tags to make it format nicely
        # if the markdown log file is displayed as HTML.
        self.info("<pre>")

        local_maxmems = [-1] * len(running_processes)
        returncodes = [None] * len(running_processes)
        proc_wrapup_text = [None] * len(running_processes)

        if not self.wait:
            self.info("</pre>")
            ids = [x.pid for x in processes]
            self.debug("Not waiting for subprocesses: " + str(ids))
            return 0, -1

        def proc_wrapup(i):
            """
            Args:
                i: internal ID number of the subprocess
            """
            returncode = processes[i].returncode
            current_pid = processes[i].pid

            info = "PID: {pid};\tCommand: {cmd};\tReturn code: {ret};\tMemory used: {mem}".format(
                pid=current_pid,
                cmd=self.running_procs[current_pid]["proc_name"],
                ret=processes[i].returncode,
                mem=display_memory(local_maxmems[i]),
            )

            # report process profile
            self._report_profile(
                self.running_procs[current_pid]["proc_name"],
                lock_file,
                time.time() - self.running_procs[current_pid]["start_time"],
                local_maxmems[i],
                current_pid,
                self.running_procs[current_pid]["args_hash"],
                self.running_procs[current_pid]["local_proc_id"],
            )

            # Remove this as a running subprocess
            self.running_procs[current_pid]["info"] = info
            self.running_procs[current_pid]["returncode"] = returncode
            self.completed_procs[current_pid] = self.running_procs[current_pid]
            del self.running_procs[current_pid]
            running_processes.remove(i)
            completed_processes.append(i)
            proc_wrapup_text[i] = info
            returncodes[i] = returncode
            return info

        sleeptime = 0.0001

        while running_processes:
            self.debug("running")
            for i in running_processes:
                local_maxmems[i] = max(local_maxmems[i], (get_mem_child_sum(processes[i])))
                self.peak_memory = max(self.peak_memory, local_maxmems[i])
                self.debug(processes[i])
                if not self._attend_process(processes[i], sleeptime):
                    proc_wrapup_text[i] = proc_wrapup(i)

            # the sleeptime is extremely short at the beginning and gets longer exponentially
            # (+ constant to prevent copious checks at the very beginning)
            # = more precise mem tracing for short processes
            sleeptime = min((sleeptime + 0.25) * 3, 60 / len(processes))

        # All jobs are done, print a final closing and job info
        info = (
            "Elapsed time: " + str(datetime.timedelta(seconds=self.time_elapsed(start_time))) + "."
        )
        info += " Running peak memory: {pipe}.".format(pipe=display_memory(self.peak_memory))
        # if len(proc_wrapup_text) == 1:
        # info += " {}".format(proc_wrapup_text[0])

        for i in completed_processes:
            info += "  \n  {}".format(self.completed_procs[processes[i].pid]["info"])

        info += "\n"  # finish out the
        self.info("</pre>")
        self.info("Command completed. {info}".format(info=info))

        for i, rc in enumerate(returncodes):
            if rc != 0:
                proc_name = self.completed_procs.get(processes[i].pid, {}).get("proc_name", cmd)
                msg = (
                    "Subprocess returned nonzero result (return code: {rc}). "
                    "Failed command: `{proc_name}`. "
                    "Check the output above for details. "
                    "To allow this step to fail without stopping the pipeline, use nofail=True.".format(
                        rc=rc, proc_name=proc_name
                    )
                )
                self._triage_error(SubprocessError(msg, returncode=rc, cmd=proc_name), nofail)

        return [returncodes, local_maxmems]

    def process_counter(self) -> int | str:
        """Increment and return the process counter for logging.

        Returns "Nf" (e.g. "3f") for commands inside a follow function,
        or the incremented integer count otherwise.
        """
        try:
            if self.in_follow:
                return str(self.proc_count) + "f"
            else:
                self.proc_count += 1
                return self.proc_count
        except AttributeError:
            self.proc_count += 1
            return self.proc_count

    ###################################
    # Waiting functions
    ###################################

    def _wait_for_process(self, p: Any, shell: bool = False) -> list[int | float]:
        """
        Debug function used in unit tests.

        Args:
            p: A subprocess.Popen process.
            shell (bool): If command requires should be run in its own shell. Optional. Default: False.
        """
        local_maxmem = -1
        sleeptime = 0.5
        while p.poll() is None:
            if not shell:
                local_maxmem = max(local_maxmem, self._memory_usage(p.pid) / 1e6)
                # print("int.maxmem (pid:" + str(p.pid) + ") " + str(local_maxmem))
            time.sleep(sleeptime)
            sleeptime = min(sleeptime + 5, 60)

        self.peak_memory = max(self.peak_memory, local_maxmem)

        del self.running_procs[p.pid]

        info = "Process " + str(p.pid) + " returned: (" + str(p.returncode) + ")."
        if not shell:
            info += " Peak memory: (Process: " + str(round(local_maxmem, 3)) + "GB;"
            info += " Pipeline: " + str(round(self.peak_memory, 3)) + "GB)\n"

        self.info(info + "\n")
        if p.returncode != 0:
            raise Exception(
                "Process {pid} returned nonzero result (return code: {rc}).".format(
                    pid=p.pid, rc=p.returncode
                )
            )
        return [p.returncode, local_maxmem]

    def _wait_for_lock(self, lock_file: str) -> None:
        """
        Just sleep until the lock_file does not exist or a lock_file-related dynamic recovery flag is spotted

        Args:
            lock_file (str): Lock file to wait upon.
        """
        sleeptime = 0.5
        first_message_flag = False
        long_message_flag = False
        dot_count = 0
        totaltime = 0
        recover_file = self._recoverfile_from_lockfile(lock_file)
        while os.path.isfile(lock_file):
            if first_message_flag is False:
                self.timestamp("Waiting for file lock: " + lock_file)
                self.warning(
                    "Lock file: {lock}. Another pipeline process may be running this step, "
                    "or a previous run crashed and left a stale lock. Options:\n"
                    "  1. Wait (current behavior) -- the lock may be released by another process.\n"
                    "  2. Restart with recover=True (CLI: -R) to override locks, discard partial output, and re-run.\n"
                    "  3. Manually delete the lock file if you are certain no other process is running:\n"
                    "     rm {lock}".format(lock=lock_file)
                )
                # self._set_status_flag(WAIT_FLAG)
                self.pipestat.set_status(
                    record_identifier=self._pipestat_manager.record_identifier,
                    status_identifier="waiting",
                )
                first_message_flag = True
            else:
                sys.stdout.write(".")
                dot_count = dot_count + 1
                if dot_count % 60 == 0:
                    self.info("")  # linefeed
            # prevents the issue of pypiper waiting for the lock file to be gone infinitely
            # in case the recovery flag is sticked by other pipeline when it's interrupted
            if os.path.isfile(recover_file):
                sys.stdout.write(" Dynamic recovery flag found")
                break
            time.sleep(sleeptime)
            totaltime += sleeptime
            sleeptime = min(sleeptime + 2.5, 60)
            if sleeptime > 3600 and not long_message_flag:
                long_message_flag = True

        if first_message_flag:
            self.timestamp("File unlocked.")
            # self._set_status_flag(RUN_FLAG)
            self.pipestat.set_status(
                record_identifier=self._pipestat_manager.record_identifier,
                status_identifier="running",
            )

    ###################################
    # Logging functions
    ###################################

    def debug(self, msg: str, *args: Any, **kwargs: Any) -> None:
        self._logger.debug(msg, *args, **kwargs)

    def info(self, msg: str, *args: Any, **kwargs: Any) -> None:
        self._logger.info(msg, *args, **kwargs)

    def warning(self, msg: str, *args: Any, **kwargs: Any) -> None:
        self._logger.warning(msg, *args, **kwargs)

    def error(self, msg: str, *args: Any, **kwargs: Any) -> None:
        self._logger.error(msg, *args, **kwargs)

    def critical(self, msg: str, *args: Any, **kwargs: Any) -> None:
        self._logger.critical(msg, *args, **kwargs)

    def fatal(self, msg: str, *args: Any, **kwargs: Any) -> None:
        self._logger.fatal(msg, *args, **kwargs)

    def timestamp(
        self,
        message: str = "",
        checkpoint: str | None = None,
        finished: bool = False,
        raise_error: bool = True,
    ) -> None:
        """Log a message with current time and elapsed time since last timestamp.

        Example:
            pm.timestamp("### Alignment")
            pm.timestamp("Reads trimmed", checkpoint="trim", finished=True)

        Messages starting with "###" are formatted as headings with surrounding
        newlines. If checkpoint is provided, creates a checkpoint file that
        enables start/stop control for pipeline reruns.

        Checkpoint semantics: with finished=False (default), the checkpoint
        marks the *start* of a stage -- used with --stop-before to halt
        before a stage runs. With finished=True, it marks *completion* of
        a stage -- used with --stop-after to halt after a stage finishes.
        A typical pattern is to call timestamp() with just a checkpoint name
        before each stage (finished=False is implied), and pypiper handles
        the rest.

        Args:
            message: Message to log with the timestamp.
            checkpoint: Stage name for checkpoint file creation.
            finished: True if the checkpoint stage just completed (retrospective);
                False (default) if the stage is just starting (prospective).
            raise_error: Whether to raise PipelineHalt when a stop point is reached.
        """

        # Halt if the manager's state has been set such that this call
        # should halt the pipeline.
        if self.halt_on_next:
            self.halt(checkpoint, finished, raise_error=raise_error)

        # Determine action to take with respect to halting if needed.
        if checkpoint:
            if finished:
                # Write the file.
                self._checkpoint(checkpoint)
                self.prev_checkpoint = checkpoint
                self.curr_checkpoint = None
            else:
                self.prev_checkpoint = self.curr_checkpoint
                self.curr_checkpoint = checkpoint
                self._checkpoint(self.prev_checkpoint)
            # Handle the two halting conditions.
            if (finished and checkpoint == self.stop_after) or (
                not finished and checkpoint == self.stop_before
            ):
                self.halt(checkpoint, finished, raise_error=raise_error)
            # Determine if we've started executing.
            elif checkpoint == self.start_point:
                self._active = True
            # If this is a prospective checkpoint, set the current checkpoint
            # accordingly and whether we should halt the pipeline on the
            # next timestamp call.
            if not finished and checkpoint == self.stop_after:
                self.halt_on_next = True

        elapsed = self.time_elapsed(self.last_timestamp)
        t = time.strftime("%m-%d %H:%M:%S")
        if checkpoint is None:
            msg = "{m} ({t}) elapsed: {delta_t} _TIME_".format(m=message, t=t, delta_t=elapsed)
        else:
            msg = "{m} ({t}) ({status} {stage}) elapsed: {delta_t} _TIME_".format(
                m=message,
                t=t,
                status="finished" if finished else "starting",
                stage=checkpoint,
                delta_t=elapsed,
            )
        if re.match("^###", message):
            msg = "\n{}\n".format(msg)
        self.info(msg)
        self.last_timestamp = time.time()

    @staticmethod
    def time_elapsed(time_since: float) -> float:
        """Return seconds elapsed since the given time (from time.time()).

        Returns:
            Elapsed seconds as float, rounded to nearest integer.
        """
        return round(time.time() - time_since, 0)

    def _report_profile(
        self,
        command: str,
        lock_name: str | None,
        elapsed_time: float,
        memory: float,
        pid: int,
        args_hash: str | None,
        proc_count: int | str,
    ) -> None:
        """
        Writes a string to self.pipeline_profile_file.
        """
        rel_lock_name = (
            lock_name if lock_name is None else os.path.relpath(lock_name, self.outfolder)
        )
        message_raw = (
            str(pid)
            + "\t"
            + str(args_hash)
            + "\t"
            + str(proc_count)
            + "\t"
            + str(datetime.timedelta(seconds=round(elapsed_time, 2)))
            + "\t "
            + str(round(memory, 4))
            + "\t"
            + str(command)
            + "\t"
            + str(rel_lock_name)
        )
        with open(self.pipeline_profile_file, "a") as myfile:
            myfile.write(message_raw + "\n")

    def report_result(
        self,
        key: str,
        value: Any,
        nolog: bool = False,
        result_formatter: Callable | None = None,
        force_overwrite: bool = True,
    ) -> Any:
        """Report a key-value result to the pipeline stats file via pipestat.

        Example:
            pm.report_result("aligned_reads", 1500000)
            pm.report_result("alignment_rate", 0.95, nolog=True)

        Args:
            key: Result name (must match schema if schema validation is on).
            value: Result value (str, int, float, dict, etc.).
            nolog: If True, suppress logging the result to the pipeline log.
            result_formatter: Custom formatter callable. Default: markdown.
            force_overwrite: Overwrite existing results. Default: True.

        Returns:
            Formatted result string(s) from pipestat.
        """
        # keep the value in memory:
        self.stats_dict[key] = value

        rf = result_formatter or self.pipestat_result_formatter

        reported_result = self.pipestat.report(
            values={key: value},
            record_identifier=self.pipestat_record_identifier,
            result_formatter=rf,
            force_overwrite=force_overwrite,
        )

        if not nolog:
            if isinstance(
                reported_result, bool
            ):  # Pipestat can return False if results are NOT reported.
                self.info("Result successfully reported? " + str(reported_result))
            else:
                for r in reported_result:
                    self.info(r)

        return reported_result

    def report_object(
        self,
        key: str,
        filename: str,
        anchor_text: str | None = None,
        anchor_image: str | None = None,
        annotation: str | None = None,
        nolog: bool = False,
        result_formatter: Callable | None = None,
        force_overwrite: bool = True,
    ) -> None:
        """Report a file/image result via pipestat.

        Deprecated: use pm.pipestat.report() directly instead:
            pm.pipestat.report(values={"peak_plot": {"path": "peaks.png"}})

        Args:
            key: Result name.
            filename: Path to the file (relative to output dir).
            anchor_text: Link text or caption. Defaults to key.
            anchor_image: Path to thumbnail image (.png/.jpg).
            annotation: Annotation string. Defaults to pipeline name.
            nolog: Suppress logging.
            result_formatter: Custom formatter callable.
            force_overwrite: Overwrite existing results.

        Returns:
            Formatted result string(s) from pipestat.
        """
        warnings.warn(
            "This function may be removed in future release. "
            "The recommended way to report pipeline results is using PipelineManager.pipestat.report().",
            category=DeprecationWarning,
        )
        rf = result_formatter or self.pipestat_result_formatter
        # Default annotation is current pipeline name.
        annotation = str(annotation or self.name)
        # In case the value is passed with trailing whitespace.
        filename = str(filename).strip()
        if anchor_text:
            anchor_text = str(anchor_text).strip()
        else:
            anchor_text = str(key).strip()
        # better to use a relative path in this file
        # convert any absolute paths into relative paths

        values = {
            "path": filename,
            "thumbnail_path": anchor_image,
            "title": anchor_text,
            "annotation": annotation,
        }
        val = {key: values}

        reported_result = self.pipestat.report(
            values=val,
            record_identifier=self.pipestat_record_identifier,
            result_formatter=rf,
            force_overwrite=force_overwrite,
        )

        if not nolog:
            if isinstance(
                reported_result, bool
            ):  # Pipestat can return False if results are NOT reported.
                self.info("Result successfully reported? " + str(reported_result))
            else:
                for r in reported_result:
                    self.info(r)

    def _safe_write_to_file(self, file: str, message: str) -> None:
        """
        Writes a string to a file safely (with file locks).
        """
        warnings.warn(
            "This function may be removed in future release. "
            "The recommended way to report pipeline results is using PipelineManager.pipestat.report().",
            category=DeprecationWarning,
        )
        target = file
        lock_name = _make_lock_name(target, self.outfolder)
        lock_file = self._make_lock_path(lock_name)

        while True:
            if os.path.isfile(lock_file):
                self._wait_for_lock(lock_file)
            else:
                try:
                    self.locks.append(lock_file)
                    self._create_file_racefree(lock_file)
                except OSError as e:
                    if e.errno == errno.EEXIST:
                        self.warning(
                            "Lock file appeared between check and creation (race condition). "
                            "Re-checking. This is normal in concurrent pipeline execution."
                        )
                        continue  # Go back to start

                # Proceed with file writing
                with open(file, "a") as myfile:
                    myfile.write(message + "\n")

                os.remove(lock_file)
                self.locks.remove(lock_file)

                # If you make it to the end of the while loop, you're done
                break

    def _report_command(self, cmd: str, procs: list[int] | None = None) -> None:
        """
        Writes a command to both stdout and to the commands log file
        (self.pipeline_commands_file).

        Args:
            cmd (str): command to report
            procs (str | list[str]): process numbers for processes in the command
        """
        if isinstance(procs, list):
            procs = ",".join(map(str, procs))
        if procs:
            line = "\n> `{cmd}` ({procs})".format(cmd=str(cmd), procs=procs)
        else:
            line = "\n> `{cmd}`".format(cmd=str(cmd))

        # Print line to stdout
        self.info(line)

        # And also add to commands file

        cmd_line = "{cmd}\n".format(cmd=str(cmd))
        with open(self.pipeline_commands_file, "a") as myfile:
            myfile.write(cmd_line)

    ###################################
    # Filepath functions
    ###################################

    @staticmethod
    def _create_file(file: str) -> None:
        """
        Creates a file, but will not fail if the file already exists.
        This is vulnerable to race conditions; use this for cases where it
        doesn't matter if this process is the one that created the file.

        Args:
            file (str): File to create.
        """
        with open(file, "w") as fout:
            fout.write("")

    @staticmethod
    def _create_file_racefree(file: str) -> None:
        """
        Creates a file, but fails if the file already exists.

        This function will thus only succeed if this process actually creates
        the file; if the file already exists, it will cause an OSError,
        solving race conditions.

        Args:
            file (str): File to create.
        """
        write_lock_flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
        fd = os.open(file, write_lock_flags)
        os.close(fd)

    @staticmethod
    def _ensure_lock_prefix(lock_name_base: str) -> str:
        """Ensure that an alleged lock file is correctly prefixed."""
        return (
            lock_name_base
            if lock_name_base.startswith(LOCK_PREFIX)
            else LOCK_PREFIX + lock_name_base
        )

    def _make_lock_path(self, lock_name_base: str) -> str:
        """
        Create path to lock file with given name as base.

        Args:
            lock_name_base (str): Lock file name, designed to not be prefixed
                with the lock file designation, but that's permitted.

        Returns:
            str: Path to the lock file.
        """

        # For lock prefix validation, separate file name from other path
        # components, as we care about the name prefix not path prefix.
        base, name = os.path.split(lock_name_base)

        lock_name = self._ensure_lock_prefix(name)
        if base:
            lock_name = os.path.join(base, lock_name)
        return _pipeline_filepath(self, filename=lock_name)

    def _recoverfile_from_lockfile(self, lockfile: str) -> str:
        """
        Create path to recovery file with given name as base.

        Args:
            lockfile (str): Name of file on which to base this path,
                perhaps already prefixed with the designation of a lock file.

        Returns:
            str: Path to recovery file.
        """
        # Require that the lock file path be absolute, or at least relative
        # and starting with the pipeline output folder.
        if not (os.path.isabs(lockfile) or lockfile.startswith(self.outfolder)):
            lockfile = self._make_lock_path(lockfile)
        return lockfile.replace(LOCK_PREFIX, "recover." + LOCK_PREFIX)

    @staticmethod
    def make_sure_path_exists(path: str) -> None:
        """Create all directories in a path, no error if they already exist.

        Args:
            path: Directory path to create.
        """
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

    ###################################
    # Pipeline stats calculation helpers
    ###################################

    def _refresh_stats(self) -> None:
        """
        Loads up the stats yaml created for this pipeline run and reads
        those stats into memory
        """

        if os.path.isfile(self.pipeline_stats_file):
            data = load_yaml(filepath=self.pipeline_stats_file)

            for key, value in data[self._pipestat_manager.pipeline_name][
                self._pipestat_manager.pipeline_type
            ][self._pipestat_manager.record_identifier].items():
                self.stats_dict[key] = value

    def get_stat(self, key: str) -> Any:
        """Retrieve a previously reported stat, loading from disk if needed.

        Example:
            trimmed = pm.get_stat("trimmed_reads")
            rate = pm.get_stat("aligned_reads") / trimmed

        Checks in-memory cache first, then reads from the stats YAML file.
        Returns None with a warning if the stat is not found.

        Args:
            key: Stat name to retrieve.

        Returns:
            The stat value, or None if not found.
        """

        try:
            return self.stats_dict[key]
        except KeyError:
            self._refresh_stats()
            try:
                return self.stats_dict[key]
            except KeyError:
                self.warning("Missing stat '{}'".format(key))
                return None

    ###################################
    # Pipeline termination functions
    ###################################

    # TODO: support checkpointing via stage, basic name, function, or filename (prohibit filepath)

    def _checkpoint(self, stage: Any) -> bool:
        """
        Decide whether to stop processing of a pipeline. This is the hook

        A pipeline can report various "checkpoints" as sort of status markers
        that designate the logical processing phase that's just been completed.
        The initiation of a pipeline can preordain one of those as a "stopping
        point" that when reached, should stop the pipeline's execution.

        Args:
            stage (pypiper.Stage | str): Pipeline processing stage/phase just completed.

        Returns:
            bool: Whether a checkpoint was created (i.e., whether it didn't
                already exist)

        Raises:
            ValueError: If the stage is specified as an absolute filepath,
                and that path indicates a location that's not immediately within
                the main output folder, raise a ValueError.
        """

        # For null stage, short-circuit and indicate no file write.
        # This handles case in which we're timestamping prospectively and
        # previously weren't in a stage.
        if stage is None:
            return False

        try:
            is_checkpoint = stage.checkpoint
        except AttributeError:
            # Maybe we have a raw function, not a stage.
            if hasattr(stage, "__call__"):
                stage = stage.__name__
            else:
                # Maybe we have a stage name not a Stage.
                # In that case, we can proceed as-is, with downstream
                # processing handling Stage vs. stage name disambiguation.
                # Here, though, warn about inputs that appear filename/path-like.
                # We can't rely on raw text being a filepath or filename,
                # because that would ruin the ability to pass stage name rather
                # than actual stage. We can issue a warning message based on the
                # improbability of a stage name containing the '.' that would
                # be expected to characterize the extension of a file name/path.
                base, ext = os.path.splitext(stage)
                if ext and "." not in base:
                    self.warning(
                        "WARNING: '{}' looks like it may be the name or path of "
                        "a file; for such a checkpoint, use touch_checkpoint.".format(stage)
                    )
        else:
            if not is_checkpoint:
                self.warning("Not a checkpoint: {}".format(stage))
                return False
            stage = stage.name

        self.info("Checkpointing: '{}'".format(stage))
        if os.path.isabs(stage):
            check_fpath = stage
        else:
            check_fpath = _checkpoint_filepath(stage, pm=self)
        return self._touch_checkpoint(check_fpath)

    def _touch_checkpoint(self, check_file: str) -> bool:
        """
        Alternative way for a pipeline to designate a checkpoint.

        Args:
            check_file (str): Name or path of file to use as checkpoint.

        Returns:
            bool: Whether a file was written (equivalent to whether the
                checkpoint file already existed).

        Raises:
            ValueError: Raise a ValueError if the argument provided as the
                checkpoint file is an absolute path and that doesn't correspond
                to a location within the main output folder.
        """
        if os.path.isabs(check_file):
            folder, _ = os.path.split(check_file)
            # For raw string comparison, ensure that each path
            # bears the final path separator.
            other_folder = os.path.join(folder, "")
            this_folder = os.path.join(self.outfolder, "")
            if other_folder != this_folder:
                errmsg = (
                    "Path provided as checkpoint file isn't in pipeline "
                    "output folder. '{}' is not in '{}'".format(check_file, self.outfolder)
                )
                raise ValueError(errmsg)
            fpath = check_file
        else:
            fpath = _pipeline_filepath(self, filename=check_file)

        # Create/update timestamp for checkpoint, but base return value on
        # whether the action was a simple update or a novel creation.
        already_exists = os.path.isfile(fpath)
        open(fpath, "w").close()
        action = "Updated" if already_exists else "Created"
        self.info("{} checkpoint file: '{}'".format(action, fpath))

        return already_exists

    def __enter__(self):
        """Support use as a context manager.

        Example:
            with PipelineManager("test", "output/") as pm:
                pm.run("echo hello", target="output/hello.txt")
            # stop_pipeline() called automatically on clean exit
            # fail_pipeline() called automatically on exception
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Clean up pipeline on context manager exit.

        Calls stop_pipeline() on clean exit, fail_pipeline() on exception.
        """
        if exc_type is None:
            self.stop_pipeline()
        else:
            try:
                self.fail_pipeline(exc_val)
            except BaseException:
                pass  # fail_pipeline re-raises; suppress to let original propagate
        return False  # Never swallow the exception

    def complete(self) -> None:
        """Mark the pipeline as successfully completed and finalize.

        Example:
            pm.complete()

        Records elapsed time and success timestamp, runs cleanup of
        intermediate files, and sets the pipeline status flag to 'completed'.
        """
        self.stop_pipeline(status=COMPLETE_FLAG)

    def fail_pipeline(self, exc: Exception, dynamic_recover: bool = False) -> None:
        """Stop the pipeline with a failure status and raise the given exception.

        Example:
            pm.fail_pipeline(Exception("BAM file is empty"))

        Terminates running subprocesses, writes cleanup script (but does
        not delete intermediate files), and sets status to 'failed'.

        Args:
            exc: Exception to raise after cleanup.
            dynamic_recover: Create recovery flag files alongside each active
                lock file. These flags signal a waiting pipeline to proceed
                rather than waiting indefinitely for a lock that will never
                be released. Used for job termination (e.g. SIGTERM from a
                cluster scheduler) rather than code errors.
        """
        self._stopped = True
        sys.stdout.flush()
        self._terminate_running_subprocesses()

        if dynamic_recover:
            # job was terminated, not failed due to a bad process.
            # flag this run as recoverable.
            if len(self.locks) < 1:
                # If there is no process locked, then recovery will be automatic.
                self.info("No locked process. Dynamic recovery will be automatic.")
            # make a copy of self.locks to iterate over since we'll be clearing them as we go
            # set a recovery flag for each lock.
            for lock_file in self.locks[:]:
                recover_file = self._recoverfile_from_lockfile(lock_file)
                self.info("Setting dynamic recover file: {}".format(recover_file))
                self._create_file(recover_file)
                self.locks.remove(lock_file)

        # Produce cleanup script
        self._cleanup(dry_run=True)

        # Finally, set the status to failed and close out with a timestamp
        if not self._failed:  # and not self._completed:
            self.timestamp("### Pipeline failed at: ")
            total_time = datetime.timedelta(seconds=self.time_elapsed(self.starttime))
            self.info("Total time: " + str(total_time))
            self.info("Failure reason: " + str(exc))
            self.pipestat.set_status(
                record_identifier=self._pipestat_manager.record_identifier,
                status_identifier="failed",
            )

        if isinstance(exc, str):
            exc = RuntimeError(exc)

        raise exc

    def halt(
        self, checkpoint: str | None = None, finished: bool = False, raise_error: bool = True
    ) -> None:
        """Pause the pipeline before its natural completion point.

        Example:
            pm.halt(checkpoint="alignment", finished=True)

        Sets status to 'paused' and raises PipelineHalt if raise_error is True.

        Args:
            checkpoint: Name of stage just reached or completed.
            finished: Whether the indicated stage was just completed.
            raise_error: Raise PipelineHalt to stop execution.
        """
        self.stop_pipeline(PAUSE_FLAG)
        self._active = False
        if raise_error:
            raise PipelineHalt(checkpoint, finished)

    def get_elapsed_time(self) -> float:
        """Calculate total pipeline runtime from the profile file.

        Example:
            total_seconds = pm.get_elapsed_time()

        Parses the profile TSV to sum unique command runtimes (deduplicating
        reruns). Falls back to wall-clock estimate if profile is unavailable.

        Returns:
            Total runtime in seconds.
        """
        if os.path.isfile(self.pipeline_profile_file):
            df = _pd.read_csv(
                self.pipeline_profile_file,
                sep="\t",
                comment="#",
                names=PROFILE_COLNAMES,
            )
            try:
                df["runtime"] = _pd.to_timedelta(df["runtime"])
            except ValueError:
                # return runtime estimate
                # this happens if old profile style is mixed with the new one
                # and the columns do not match
                return self.time_elapsed(self.starttime)
            unique_df = df[~df.duplicated("cid", keep="last").values]
            return sum(unique_df["runtime"].apply(lambda x: x.total_seconds()))
        return self.time_elapsed(self.starttime)

    def stop_pipeline(self, status: str = COMPLETE_FLAG) -> None:
        """Terminate the pipeline, recording time/memory stats and running cleanup.

        Example:
            pm.stop_pipeline()  # defaults to 'completed' status

        This is the underlying shutdown function. It reports Time and Success
        results, removes lock files, deletes intermediate files registered
        via clean_add(), and sets the status flag. Prefer complete() for
        success or fail_pipeline() for failure -- they call this internally.

        Args:
            status: Status flag string. Default: 'completed'.
        """
        self._stopped = True
        self.pipestat.set_status(
            record_identifier=self._pipestat_manager.record_identifier,
            status_identifier=status,
        )
        self._cleanup()
        elapsed_time_this_run = str(datetime.timedelta(seconds=self.time_elapsed(self.starttime)))
        self.report_result("Time", elapsed_time_this_run, nolog=True)
        self.report_result("Success", time.strftime("%m-%d-%H:%M:%S"), nolog=True)

        self.info("\n### Pipeline completed. Epilogue")
        # print("* " + "Total elapsed time".rjust(20) + ":  "
        #       + str(datetime.timedelta(seconds=self.time_elapsed(self.starttime))))
        self.info("* " + "Elapsed time (this run)".rjust(30) + ":  " + elapsed_time_this_run)
        self.info(
            "* "
            + "Total elapsed time (all runs)".rjust(30)
            + ":  "
            + str(datetime.timedelta(seconds=round(self.get_elapsed_time())))
        )
        self.info(
            "* "
            + "Peak memory (this run)".rjust(30)
            + ":  "
            + str(round(self.peak_memory, 4))
            + " GB"
        )
        # self.info("* " + "Total peak memory (all runs)".rjust(30) + ":  " +
        #     str(round(self.peak_memory, 4)) + " GB")
        if self.halted:
            return

        t = time.strftime("%Y-%m-%d %H:%M:%S")
        self.info("* " + "Pipeline completed time".rjust(30) + ": " + t)

    def _signal_term_handler(self, signal: int, frame: Any) -> None:
        """
        TERM signal handler function: this function is run if the process
        receives a termination signal (TERM). This may be invoked, for example,
        by SLURM if the job exceeds its memory or time limits. It will simply
        record a message in the log file, stating that the process was
        terminated, and then gracefully fail the pipeline. This is necessary to
        1. set the status flag and 2. provide a meaningful error message in the
        tee'd output; if you do not handle this, then the tee process will be
        terminated before the TERM error message, leading to a confusing log
        file.
        """
        signal_type = "SIGTERM"
        self._generic_signal_handler(signal_type)

    def _generic_signal_handler(self, signal_type: str) -> None:
        """
        Function for handling both SIGTERM and SIGINT
        """
        self.info("</pre>")
        message = "Got " + signal_type + ". Failing gracefully..."
        self.timestamp(message)
        self.fail_pipeline(KeyboardInterrupt(signal_type), dynamic_recover=True)
        sys.exit(1)

        # I used to write to the logfile directly, because the interrupts
        # would first destroy the tee subprocess, so anything printed
        # after an interrupt would only go to screen and not get put in
        # the log, which is bad for cluster processing. But I later
        # figured out how to sequester the kill signals so they didn't get
        # passed directly to the tee subprocess, so I could handle that on
        # my own; hence, now I believe I no longer need to do this. I'm
        # leaving this code here as a relic in case something comes up.
        # with open(self.pipeline_log_file, "a") as myfile:
        #   myfile.write(message + "\n")

    def _signal_int_handler(self, signal: int, frame: Any) -> None:
        """
        For catching interrupt (Ctrl +C) signals. Fails gracefully.
        """
        signal_type = "SIGINT"
        self._generic_signal_handler(signal_type)

    def _exit_handler(self) -> None:
        """
        This function I register with atexit to run whenever the script is completing.
        A catch-all for uncaught exceptions, setting status flag file to failed.
        """

        # Make the cleanup file executable if it exists
        if os.path.isfile(self.cleanup_file):
            # Make the cleanup file self destruct.
            with open(self.cleanup_file, "a") as myfile:
                myfile.write("rm " + self.cleanup_file + "\n")
            os.chmod(self.cleanup_file, 0o755)

        # If the pipeline hasn't completed successfully, or already been marked
        # as failed, then mark it as failed now. During interpreter shutdown
        # (e.g., pytest cleanup), logging may fail. See https://bugs.python.org/issue11380
        if not self._has_exit_status:
            # Disable logging to avoid "I/O operation on closed file" errors
            import logging

            logging.disable(logging.CRITICAL)
            try:
                self.fail_pipeline(
                    Exception(
                        "Pipeline exited without calling stop_pipeline() or complete(). "
                        "This usually means an unhandled exception occurred earlier in the pipeline. "
                        "Check the log output above for the original error."
                    )
                )
            except Exception:
                # Silently ignore any errors during shutdown
                pass
            finally:
                logging.disable(logging.NOTSET)

        if self.tee:
            self.tee.kill()

    def _terminate_running_subprocesses(self) -> None:
        # make a copy of the list to iterate over since we'll be removing items
        for pid in self.running_procs.copy():
            proc_dict = self.running_procs[pid]

            # Close the preformat tag that we opened when the process was spawned.
            # record profile of any running processes before killing
            elapsed_time = time.time() - self.running_procs[pid]["start_time"]
            process_peak_mem = self._memory_usage(pid, container=proc_dict["container"]) / 1e6
            self._report_profile(
                self.running_procs[pid]["proc_name"],
                None,
                elapsed_time,
                process_peak_mem,
                pid,
                self.running_procs[pid]["args_hash"],
                self.running_procs[pid]["local_proc_id"],
            )
            self._kill_child_process(pid, proc_dict["proc_name"])
            del self.running_procs[pid]

    def _kill_child_process(self, child_pid: int | None, proc_name: str | None = None) -> None:
        """
        Pypiper spawns subprocesses. We need to kill them to exit gracefully,
        in the event of a pipeline termination or interrupt signal.
        By default, child processes are not automatically killed when python
        terminates, so Pypiper must clean these up manually.
        Given a process ID, this function just kills it.

        Args:
            child_pid (int): Child process id.
        """

        # When we kill process, it turns into a zombie, and we have to reap it.
        # So we can't just kill it and then let it go; we call wait

        def pskill(proc_pid, sig=signal.SIGINT):
            parent_process = psutil.Process(proc_pid)
            for child_proc in parent_process.children(recursive=True):
                child_proc.send_signal(sig)
            parent_process.send_signal(sig)

        if child_pid is None:
            return

        if proc_name:
            proc_string = " ({proc_name})".format(proc_name=proc_name)

        # First a gentle kill
        sys.stdout.flush()
        still_running = self._attend_process(psutil.Process(child_pid), 0)
        sleeptime = 0.25
        time_waiting = 0

        while still_running and time_waiting < 3:
            try:
                if time_waiting > 2:
                    pskill(child_pid, signal.SIGKILL)
                    # print("pskill("+str(child_pid)+", signal.SIGKILL)")
                elif time_waiting > 1:
                    pskill(child_pid, signal.SIGTERM)
                    # print("pskill("+str(child_pid)+", signal.SIGTERM)")
                else:
                    pskill(child_pid, signal.SIGINT)
                    # print("pskill("+str(child_pid)+", signal.SIGINT)")

            except OSError:
                # This would happen if the child process ended between the check
                # and the next kill step
                still_running = False
                time_waiting = time_waiting + sleeptime

            # Now see if it's still running
            time_waiting = time_waiting + sleeptime
            if not self._attend_process(psutil.Process(child_pid), sleeptime):
                still_running = False

        if still_running:
            # still running!?
            self.warning(
                "Child process {child_pid}{proc_string} did not terminate after repeated signals (SIGINT, SIGTERM, SIGKILL). "
                "The process may be in an uninterruptible state (e.g., waiting on I/O). "
                "You may need to manually kill it: kill -9 {child_pid}".format(
                    child_pid=child_pid, proc_string=proc_string
                )
            )
        else:
            if time_waiting > 0:
                note = "terminated after {time} sec".format(time=int(time_waiting))
            else:
                note = "was already terminated"

            msg = "Child process {child_pid}{proc_string} {note}.".format(
                child_pid=child_pid, proc_string=proc_string, note=note
            )
            self.info(msg)

    @staticmethod
    def _atexit_register(*args: Any) -> None:
        """Convenience alias to register exit functions without having to import atexit in the pipeline."""
        atexit.register(*args)

    def get_container(self, image: str, mounts: str | list[str]) -> None:
        """Start a Docker container for running pipeline commands.

        Example:
            pm.get_container("nsheff/refgenie", ["/data", "/ref"])

        Args:
            image: Docker image name (e.g. "nsheff/refgenie").
            mounts: Path or list of paths to mount into the container.
        """
        if isinstance(mounts, str):
            mounts = [mounts]
        cmd = "docker run -itd"
        for mnt in mounts:
            absmnt = os.path.abspath(mnt)
            cmd += " -v " + absmnt + ":" + absmnt
        cmd += " -v {cwd}:{cwd} --workdir={cwd}".format(cwd=os.getcwd())
        cmd += " --user={uid}".format(uid=os.getuid())
        cmd += " --volume=/etc/group:/etc/group:ro"
        cmd += " --volume=/etc/passwd:/etc/passwd:ro"
        cmd += " --volume=/etc/shadow:/etc/shadow:ro"
        cmd += " --volume=/etc/sudoers.d:/etc/sudoers.d:ro"
        cmd += " --volume=/tmp/.X11-unix:/tmp/.X11-unix:rw"
        cmd += " " + image
        container = self.checkprint(cmd).rstrip()
        self.container = container
        self.info("Using docker container: " + container)
        self._atexit_register(self.remove_container, container)

    def remove_container(self, container: str) -> None:
        """Remove a Docker container by ID.

        Args:
            container: Docker container ID string.
        """
        self.info("Removing docker container. . .")
        cmd = "docker rm -f " + container
        self.callprint(cmd)

    def clean_add(
        self, regex: str | None, conditional: bool = False, manual: bool = False
    ) -> None:
        """Register files for automatic deletion when the pipeline succeeds.

        Example:
            pm.clean_add("intermediate.bam")
            pm.clean_add("temp_*.txt", conditional=True)
            pm.clean_add(None)  # no-op, safe for unset variables

        Passing None is a no-op, allowing safe calls on variables that may
        not have been assigned in conditional branches.

        Args:
            regex: Unix glob pattern or filename to delete. None is a no-op.
            conditional: Only delete if no other pipelines are running
                (checks for absence of flag files from other pipelines).
            manual: Only add to the manual cleanup script, never auto-delete.
                Note: if the PipelineManager was created with dirty=True,
                all files are forced to manual cleanup regardless.
        """
        if regex is None:
            return

        # TODO: print this message (and several below) in debug
        # print("Adding regex to cleanup: {}".format(regex))
        if self.dirty:
            # Override the user-provided option and force manual cleanup.
            manual = True

        if not self.clean_initialized:
            # Make cleanup files relative to the cleanup script in case the result folder moves.
            with open(self.cleanup_file, "a") as myfile:
                clean_init = 'DIR="$(cd -P -- "$(dirname -- "$0")" && pwd -P)"'
                myfile.write(clean_init + "\n")
                myfile.write("cd ${DIR}\n")
                self.clean_initialized = True

        if manual:
            filenames = glob.glob(regex)
            if not filenames:
                self.info("No files match cleanup pattern: {}".format(regex))
            for filename in filenames:
                try:
                    with open(self.cleanup_file, "a") as myfile:
                        if os.path.isabs(filename):
                            relative_filename = os.path.relpath(filename, self.outfolder)
                            absolute_filename = filename
                        else:
                            relative_filename = os.path.relpath(filename, self.outfolder)
                            absolute_filename = os.path.abspath(
                                os.path.join(self.outfolder, relative_filename)
                            )
                        if os.path.isfile(absolute_filename):
                            # print("Adding file to cleanup: {}".format(filename))
                            myfile.write("rm " + relative_filename + "\n")
                        elif os.path.isdir(absolute_filename):
                            # print("Adding directory to cleanup: {}".format(filename))
                            # first, add all filenames in the directory
                            myfile.write("rm " + relative_filename + "/*\n")
                            # and the directory itself
                            myfile.write("rmdir " + relative_filename + "\n")
                        else:
                            self.info("File not added to cleanup: {}".format(relative_filename))
                except Exception as e:
                    self.error("Error in clean_add on path {}: {}".format(filename, str(e)))
        elif conditional:
            self.cleanup_list_conditional.append(regex)
        else:
            self.cleanup_list.append(regex)
            # TODO: what's the "absolute" list?
            # Remove it from the conditional list if added to the absolute list
            while regex in self.cleanup_list_conditional:
                self.cleanup_list_conditional.remove(regex)

    def _cleanup(self, dry_run: bool = False) -> None:
        """
        Cleans up (removes) intermediate files.

        You can register intermediate files, which will be deleted automatically
        when the pipeline completes. This function deletes them,
        either absolutely or conditionally. It is run automatically when the
        pipeline succeeds, so you shouldn't need to call it from a pipeline.

        Args:
            dry_run (bool): Set to True if you want to build a cleanup
                script, but not actually delete the files.
        """

        n_to_clean = len(self.cleanup_list)
        n_to_clean_cond = len(self.cleanup_list_conditional)

        if n_to_clean + n_to_clean_cond > 0:
            self.info(
                "Starting cleanup: {} files; {} conditional files for cleanup".format(
                    n_to_clean, n_to_clean_cond
                )
            )
        else:
            self.debug("No files to clean.")

        if dry_run:
            # Move all unconditional cleans into the conditional list
            if n_to_clean > 0:
                combined_list = self.cleanup_list_conditional + self.cleanup_list
                self.cleanup_list_conditional = combined_list
                self.cleanup_list = []

        if n_to_clean > 0:
            self.info("\nCleaning up flagged intermediate files. . .")
            for expr in self.cleanup_list:
                self.debug("Removing glob: " + expr)
                try:
                    # Expand regular expression
                    files = glob.glob(expr)
                    # Remove entry from cleanup list
                    while files in self.cleanup_list:
                        self.cleanup_list.remove(files)
                    # and delete the files
                    for file in files:
                        if os.path.isfile(file):
                            self.debug("`rm {}`".format(file))
                            os.remove(os.path.join(file))
                        elif os.path.isdir(file):
                            self.debug("`rmdir {}`".format(file))
                            os.rmdir(os.path.join(file))
                except Exception:
                    pass

        if n_to_clean_cond > 0:
            run_flag = _flag_name(RUN_FLAG)
            flag_files = [
                fn
                for fn in glob.glob(self.outfolder + _flag_name("*"))
                if COMPLETE_FLAG not in os.path.basename(fn)
                and not "{}_{}_{}".format(
                    self._pipestat_manager.pipeline_name,
                    self.pipestat_record_identifier,
                    run_flag,
                )
                == os.path.basename(fn)
            ]
            if len(flag_files) == 0 and not dry_run:
                self.info("\nCleaning up conditional list. . .")
                for expr in self.cleanup_list_conditional:
                    self.debug("\nRemoving glob: " + expr)
                    try:
                        files = glob.glob(expr)
                        while files in self.cleanup_list_conditional:
                            self.cleanup_list_conditional.remove(files)
                        for file in files:
                            if os.path.isfile(file):
                                self.debug("`rm {}`".format(file))
                                os.remove(os.path.join(file))
                            elif os.path.isdir(file):
                                self.debug("`rmdir {}`".format(file))
                                os.rmdir(os.path.join(file))
                    except Exception:
                        pass
            else:
                self.info(
                    "\nConditional flag found: " + str([os.path.basename(i) for i in flag_files])
                )
                self.info(
                    "\nThese conditional files were left in place:\n\n- "
                    + "\n- ".join(self.cleanup_list_conditional)
                )
                # Produce a cleanup script.
                no_cleanup_script = []
                for cleandir in self.cleanup_list_conditional:
                    try:
                        items_to_clean = glob.glob(cleandir)
                        for clean_item in items_to_clean:
                            with open(self.cleanup_file, "a") as clean_script:
                                if os.path.isfile(clean_item):
                                    clean_script.write("rm " + clean_item + "\n")
                                elif os.path.isdir(clean_item):
                                    clean_script.write("rmdir " + clean_item + "\n")
                    except Exception:
                        no_cleanup_script.append(cleandir)
                if no_cleanup_script:
                    self.warning(
                        "\n\nCould not produce cleanup script for item(s):\n\n- "
                        + "\n- ".join(no_cleanup_script)
                    )

    def _memory_usage(
        self, pid: int | str = "self", category: str = "hwm", container: str | None = None
    ) -> float | int:
        """
        Memory usage of the process in kilobytes.

        Args:
            pid (str): Process ID of process to check
            category (str): Memory type to check. 'hwm' for high water mark.
        """
        if container:
            # TODO: Put some debug output here with switch to Logger
            # since this is relatively untested.
            cmd = "docker stats " + container + " --format '{{.MemUsage}}' --no-stream"
            mem_use_str = subprocess.check_output(cmd, shell=True).decode()

            mem_num = re.findall(r"[\d\.]+", mem_use_str.split("/")[0])[0]
            mem_scale = re.findall(r"[A-Za-z]+", mem_use_str.split("/")[0])[0]

            mem_num = float(mem_num)
            if mem_scale == "GiB":
                return mem_num * 1e6
            elif mem_scale == "MiB":
                return mem_num * 1e3
            elif mem_scale == "KiB":
                return mem_num
            else:
                # What type is this?
                return 0

        # Thanks Martin Geisler:
        status = None
        result = {"peak": 0, "rss": 0, "hwm": 0}

        try:
            # This will only work on systems with a /proc file system
            # (like Linux).
            # status = open('/proc/self/status')
            proc_spot = "/proc/%s/status" % pid
            status = open(proc_spot)
            for line in status:
                parts = line.split()
                key = parts[0][2:-1].lower()
                if key in result:
                    result[key] = int(parts[1])
        except Exception:
            return 0

        finally:
            if status is not None:
                status.close()
        # print(result[category])
        return result[category]

    def _triage_error(self, e: Exception, nofail: bool) -> None:
        """Print a message and decide what to do about an error."""
        if not nofail:
            self.fail_pipeline(e)
        elif self._failed:
            self.info(
                "This step has nofail=True, but the pipeline was already marked as failed "
                "by a previous error, so this error is also raised. "
                "Fix the earlier failure first."
            )
            raise e
        else:
            self.error(e)
            self.error(
                "Subprocess returned nonzero result, but pipeline is continuing because nofail=True. "
                "To make this a fatal error, remove the nofail=True argument from this run() call."
            )
        # TODO: return nonzero, or something. . .?
