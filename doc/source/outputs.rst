
Outputs
=========================

Assume you name your pipeline `PIPE` (by `passing name="PIPE"` to the PipelineManager constructor), by default, your ``PipelineManager`` will produce the following outputs automatically (in additional to any output created by the actual pipeline commands you run):

* **PIPE_log.md** - the log starts with a bunch of useful information about your run: a starting timestamp, version numbers of the pipeline and pypiper, a declaration of all arguments passed to the pipeline, the compute host, etc. Then, all output sent to screen is automatically logged to this file, providing a complete record of your run without requiring you to write **any** logging code on your own.

* **PIPE_status.flag** - As the pipeline runs, it produces a flag in the output directory, which can be either `PIPE_running.flag`, `PIPE_failed.flag`, or `PIPE_completed.flag`. These flags make it easy to assess the current state of running pipelines for individual samples, and for many samples in a project simultaneously (the `flagCheck.sh` script produces a summary of all pipeline runs in an output folder).

* **PIPE_stats.md** - any results reported using ``report_result()`` are saved as key-value pairs in this file, for easy parsing (we also have a post-pipeline script to aggregate these results into a summary table).

* **PIPE_profile.md** - A profile log file that provides, for every process run by the pipeline, 3 items: 1) the process name; 2) the clock time taken by the process; and 3) the high water mark memory used by the process. This file makes it easy to profile pipelines for memory and time resources.

* **PIPE_commands.md** - Pypiper produces a log file containing all the commands run by the pipeline, verbatim. These are also included in the main log.

Multiple pipelines can easily be run on the same sample, using the same output folder (and possibly sharing intermediate files), as the result outputs will be identifiable by the `PIPE` identifier.

These files are `markdown <https://daringfireball.net/projects/markdown/>`_ making it easy to read either in text format, or to quickly convert to a pretty format like HTML.
