# Pypiper
A lightweight python toolkit for gluing together restartable, robust command line pipelines.

The public-facing permalink for the project is http://databio.com/pypiper. The source code, tutorials, are hosted at http://github.com/epigen/pypiper.

# Introduction

Pypiper helps you produce pipelines. 

The target user of Pypiper is a computational scientist comfortable on the command line, who has something like a `bash` script that would benefit from a layer of "handling code". Pypiper helps you convert that set of shell commands into a production-scale workflow, automatically handling the annoying details to make your pipeline robust and restartable, with minimal learning curve.

Pypiper *does not* handle any sort of cluster job submission, resource requesting, or parallel dependency management. You can use your current setup for those things, and use Pypiper just to produce a robust, restartable, and logged procedural pipeline.

# Benefits of using Pypiper

* Restartability - Commands check for their targets and only run if the target needs to be created, much like a `makefile`, making the pipeline pick up where it left off in case it needs to be restarted or extended.
* Pipeline integrity protection - PyPiper uses file locking to ensure that multiple pipeline runs will not interfere with one another -- even if the steps are identical and produce the same files. One run will seemless wait for the other, making it possible to share steps seemlessly across pipelines.
* Memory use monitoring - Processes are polled for high water mark memory use, allowing you to more accurately guage your future memory requirements.
* Easy job status monitoring - Pypiper uses status flag files to make it possible to assess the current state (`running`, `failed`, or `completed`) of hundreds of jobs simultaneously.
* Robust error handling - Pypiper closes pipelines gracefully on interrupt or termination signals, converting the status to `failed`. By default, a process that returns a nonzero value halts the pipeline, unlike in bash, where by default the pipeline would continue using an incomplete or failed result. This behavior can be overridden as desired with a single parameter.
* Logging - Pypiper automatically records the output of your pipeline and its subprocesses, and provides copious information on pipeline initiation, as well as easy timestamping.
* Easy result reports - Pypiper provides functions to put key-value pairs into an easy-to-parse stats file, making it easy to summarize your pipeline results.
* Simplicity - It should only take you 15 minutes to run your first pipeline. The basic documentation is just a few pages long. The codebase itself is also only a few thousand lines of code, making it very lightweight.


# Your first pipeline

Using Pypiper is simple. First, import pypiper, and then create a new Pypiper object:

```{python}
import pypiper, os
outfolder = "pipeline_output/" # Choose a folder for your results
pipeline = pypiper.Pypiper(name="my_pipeline", outfolder=outfolder)
```

Now, the workhorse of Pypiper is the `call_lock()` function. Essentially, you just create a shell command as a string in python, and then pass it to `call_lock()`. 

```
target = os.path.join(outfolder, "outfile.txt")
command cmd = "shuf -i 1-500000000 -n 10000000 > " + target
pipeline.call_lock(command, target, shell=True)
```

There's more information about `shell=True` processes below. Now string together whatever commands your pipeline requires, and then at the end, terminate the pipeline so it gets flagged as successfully completed:

```
pipeline.stop_pipeline()
```

That's it! By running this command through `call_lock()` instead of directly in bash, you get a robust, logged, restartable pipeline manager for free!

# Tutorials

To see an example of a simple pipline, look in the `example_pipelines` folder in this respository, which are thoroughly commented to act as vignettes. This is the best way to learn how to use Pyipiper.

* [example_pipelines/basic.py](example_pipelines/basic.py) - How to construct a command and use `call_lock()`.
* `example_pipelines/advanced.py` (under construction) - A tutorial demonstrating some more advanced features.
* `example_pipelines/bioinformatics.py` (under construction) - A tutorial showing some bioinformatics use cases.

# Outputs

Assume you name your pipeline `PIPE` (by passing name="PIPE" to the Pypiper constructor), by default, Pypiper will produce the following outputs automatically (in additional to any output created by the actual pipeline commands you run):

* `PIPE_log.md` - all output sent to screen is automatically logged to this file, including versions, start and stop time, compute host, and other details.
* `PIPE_stats.md` - any results reported using `report_result()` are saved as key-value pairs in this file, for easy parsing (we also have a post-pipeline script to aggregate these results into a summary table).
* `PIPE_status.flag` - As the pipeline runs, it produces a flag in the output directory, which can be either `PIPE_running.flag`, `PIPE_failed.flag`, or `PIPE_completed.flag`. These flags make it easy to assess the current state of running pipelines for individual samples, and for many samples in a project simultaneously (the `flagCheck.sh` script produces a summary of all pipeline runs in an output folder).
* `PIPE_profile.md` - A profile log file that provides, for every process run by the pipeline, 3 items: 1) the process name; 2) the clock time taken by the process; and 3) the high water mark memory used by the process. This file makes it easy to profile pipelines for memory and time resources.
* `PIPE_commands.md` - Pypiper produces a log file containing all the commands run by the pipeline, verbatim. These are also included in the main log.

Multiple pipelines can easily be run on the same sample, using the same output folder (and possibly sharing intermediate files), as the result outputs will be identifiable by the `PIPE` identifier.

# Toolkits

An optional feature of pypiper is the accompanying toolkits, such as the next-gen sequencing toolkit, [ngstk](pypiper/ngstk.py), which simply provides some convenient helper functions to create common commands, like converting from file formats (_e.g._ bam to fastq), merging files (_e.g._ merge_bams), counting reads, etc. These make it faster to design bioinformatics pipelines in Pypiper, but are entirely optional. Contributions of additional toolkits or functions in an existing toolkit are welcome.

# Technical Documentation
You can use `make` to generate the PyPiper documentation. Just change your working directory to `doc` and run `make` to see available documentation formats *e.g.*: `make html`. The documentation will be produced under `doc/build`.

# Testing

You can test pypiper by running `python test_pypiper.py`, which has some unit tests.

# Python process types: Shell vs direct

Since Pypiper runs all your commands from within python (using the `subprocess` python module), you need to be aware of the two types of processes that `subprocess` can handle: Shell process and direct processes. Choosing between the two is as simple as specifying `shell=True` for shell process (the default is a direct subprocess) to `call_lock()`, which is directly passed on as the `shell` parameter to `subprocess.Popen()`. The difference is that with a direct subprocess, Python just runs the process as a subprocess directly, while with a shell process, it spawns a shell process, and then uses the shell process to run your subprocess. When should you use shell processes?

**Direct process**: For most use cases, you should simply use a direct subprocess (don't pass `shell=True`) -- this has the advantage of enabling Python to monitor the memory use of the subprocess, because Python retains control over it. This is considered by the community to be the preferable way of running subprocesses in Python.

**Shell process**: You must use a `shell=True` process if you are using shell operators in your processes command. For instance, if you use an asterisk (`*`) for wildcard expansion, or a bracket (`>`) for output redirection, or a pipe (`|`) to link processes -- these are commands understood by a shell like Bash, and thus, cannot be run as direct subprocesses by Python, but must be run in a shell. `Subprocess` (and therefore `Pypiper`) can understand these commands just fine, but you lose the ability to monitor memory high water mark because Python does not have direct control over subprocesses run inside a subshell.

# Motivation
As I began to put together production-scale pipelines, I found a lot of relevant pipelining systems, but was universally disappointed. For my needs, they were all overly complex. I wanted something simple enough to quickly write and maintain a pipeline without having to learn a lot of new functions and conventions, but robust enough to handle requirements like restartability and memory usage monitoring. Everything related was either a pre-packaged pipeline for a defined purpose, or a heavy-duty development environment that was overkill for my needs. Both of these seemed to be targeted toward less experienced developers who sought structure, and neither fit my needs: I had a set of commands already in mind; I just needed a wrapper that could take that code and make it automatically restartable, logged, robust to crashing, easy to debug, and so forth.

If you need a full-blown environment that can do everything, look elsewhere. Pypiper's strength is its simplicity. If all you want is a shell-like script, but now with the power of python, and restartability, then Pypiper is for you.

### Other software

If Pypiper is not for you, check out some of these other pipelining frameworks:

* bpipe
* Anduril
* Ruffus


