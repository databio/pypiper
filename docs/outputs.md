
# Outputs explained

Assume you are using a pypiper pipeline named `PIPE` ( it passes `name="PIPE"` to the PipelineManager constructor). By default, your `PipelineManager` will produce the following outputs automatically (in addition to any output created by the actual pipeline commands you run):

* **PIPE_log.md**
	The log starts with a bunch of useful information about your run: a starting timestamp, version numbers of the pipeline and pypiper, a declaration of all arguments passed to the pipeline, the compute host, etc. Then, all output sent to screen is automatically logged to this file, providing a complete record of your run.

* **PIPE_status.flag**
	As the pipeline runs, it produces a flag in the output directory, which can be either `PIPE_running.flag`, `PIPE_failed.flag`, or `PIPE_completed.flag`. These flags make it easy to assess the current state of running pipelines for individual samples, and for many samples in a project simultaneously.

* **stats.yaml**
	Any results reported by the pipeline are saved as key-value pairs in this file, for easy parsing.

* **PIPE_profile.md**
	A profile log file that provides, for every process run by the pipeline, 3 items: 1) the process name; 2) the clock time taken by the process; and 3) the memory high water mark used by the process. This file makes it easy to profile pipelines for memory and time resources.

* **PIPE_commands.md**
	Pypiper produces a log file containing all the commands run by the pipeline, verbatim. These are also included in the main log.

Multiple pipelines can easily be run on the same sample, using the same output folder (and possibly sharing intermediate files), as the result outputs will be identifiable by the `PIPE_` identifier.

These files are [markdown](https://daringfireball.net/projects/markdown/) making it easy to read either in text format, or to quickly convert to a pretty format like HTML.
