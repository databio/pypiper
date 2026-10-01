# Your first basic pipeline

## The basic functions

Pypiper is simple but powerful. Your pipeline is a python script, let's call it `pipeline.py`. You really only need to know about 3 commands to get started:

* `PipelineManager()`: PipelineManager constructor, initializes the pipeline.
* `PipelineManager.run()`: The primary workhorse function; runs a command.
* `PipelineManager.stop_pipeline()`: Terminate the pipeline.

That's all you need to create a powerful pipeline. You can find in-depth reference documentation for each method in the API. In particular, most of Pypiper's power comes from the `run` method, which has a series of options outlined in [dedicated documentation on the run method](../advanced-run-method.md).

To write your first basic pipeline, first `import pypiper`, then specify an output folder and create a new `PipelineManager` object:

	#!/usr/bin/env python
	import pypiper, os
	outfolder = "pipeline_output/"  # Choose a folder for your results
	pm = pypiper.PipelineManager(name="my_pipeline", outfolder=outfolder)

Creating the `PipelineManger` object creates your `outfolder` and places a flag called `my_pipeline_running.flag` in the folder. It also initializes the log file (``my_pipeline_log.md``) with statistics such as time of starting, compute node, software versions, command-line parameters, etc.

Now, the workhorse of `PipelineManager` is the `run()` function. Essentially, you create a shell command as a string in python, and then pass it and its target (a file it creates) to `run()`. The target is the final output file created by your command. Let's use the built-in `shuf` command to create some random numbers and put them in a file called `outfile.txt`:

	# our command will produce this output file
	target = os.path.join(outfolder, "outfile.txt")
	command = "shuf -i 1-500000000 -n 10000000 > " + target
	pm.run(command, target)


The command (`command`) is the only required argument to `run()`. You can leave `target` empty (pass `None`). If you **do** specify a target, the command will only be run if the target file does not already exist. If you **do not** specify a target, the command will be run every time the pipeline is run. 

Now string together whatever commands your pipeline requires! At the end, terminate the pipeline so it gets flagged as successfully completed:

	pm.stop_pipeline()

That's it! By running commands through `run()` instead of directly in bash, you get a robust, logged, restartable pipeline manager for free!



## Advanced functions
Once you've mastered the basics, add in a few other commands that make debugging and development easier:

* `timestamp()`: Print message, time, and time elapsed, perhaps creating checkpoint.
* `clean_add()`: Add files (or regexs) to a cleanup list, to delete when this pipeline completes successfully.

## Reporting results 
`Pypiper` also has a really nice system for collecting and aggregating result statistics and pipeline outputs. 

* `report_result()`: Writes a line in the `stats.tsv` file with a key-value pair of some pipeline result.
* `report_object()`: Writes a line in the `objects.tsv` file with a key-value pair pointing to any file (such as an image, HTML report, or other) created by the pipeline.
* `get_stat()`:	Returns a stat that was previously reported.

## Putting them together in `basic.py`
Now, download [basic.py](https://github.com/databio/pypiper/blob/master/example_pipelines/basic.py) and run it with `python basic.py` (or, better yet, make it executable (`chmod 755 basic.py`) and then run it directly with `./basic.py`). This example is a documented vignette; so just read it and run it to get an idea of how things work. Here are the contents of `basic.py`:


```python
# %load ../example_pipelines/basic.py
#!/usr/bin/env python

"""Getting Started: A simple sample pipeline built using pypiper."""

# This is a runnable example. You can run it to see what the output
# looks like.

# First, make sure you can import the pypiper package

import os
import pypiper

# Create a PipelineManager instance (don't forget to name it!)
# This starts the pipeline.

pm = pypiper.PipelineManager(name="BASIC",
    outfolder="pipeline_output/")

# Now just build shell command strings, and use the run function
# to execute them in order. run needs 2 things: a command, and the
# target file you are creating.

# First, generate some random data

# specify target file:
tgt = "pipeline_output/test.out"

# build the command
cmd = "shuf -i 1-500000000 -n 10000000 > " + tgt

# and run with run().
pm.run(cmd, target=tgt)

# Now copy the data into a new file.
# first specify target file and build command:
tgt = "pipeline_output/copied.out"
cmd = "cp pipeline_output/test.out " + tgt
pm.run(cmd, target=tgt)

# You can also string multiple commands together, which will execute
# in order as a group to create the final target.
cmd1 = "sleep 5"
cmd2 = "touch pipeline_output/touched.out"
pm.run([cmd1, cmd2], target="pipeline_output/touched.out")

# A command without a target will run every time.
# Find the biggest line
cmd = "awk 'n < $0 {n=$0} END{print n}' pipeline_output/test.out"
pm.run(cmd, "lock.max")

# Use checkprint() to get the results of a command, and then use
# report_result() to print and log key-value pairs in the stats file:
last_entry = pm.checkprint("tail -n 1 pipeline_output/copied.out")
pm.report_result("last_entry", last_entry)


# Now, stop the pipeline to complete gracefully.
pm.stop_pipeline()

# Observe your outputs in the pipeline_output folder 
# to see what you've created.

```

# Building a pipeline to count number of reads

Here we have a more advanced bioinformatics pipeline that adds some new concepts. This is a simple script that takes an input file and returns the file size and the number of sequencing reads in that file. This example uses a function from from the built-in :doc:`NGSTk toolkit <ngstk>`. In particular, this toolkit contains a few handy functions that make it easy for a pipeline to accept inputs of various types. So, this pipeline can count the number of reads from files in ``BAM`` format, or ``fastq`` format, or ``fastq.gz`` format. You can also use the same functions from NGSTk to develop a pipeline to do more complicated things, and handle input of any of these types.

First, grab this pipeline. Download [count_reads.py](https://github.com/databio/pypiper/blob/master/example_pipelines/count_reads.py), make it executable (`chmod 755 count_reads.py`), and then run it with `./count_reads.py`). 

You can grab a few small data files in the [microtest repository](https://github.com/epigen/microtest/tree/master/config). Run a few of these files like this:

```
./count_reads.py -I ~/code/microtest/data/rrbs_PE_R1.fastq.gz -O $HOME -S sample1
./count_reads.py -I ~/code/microtest/data/rrbs_PE_fq_R1.fastq -O $HOME -S sample2
./count_reads.py -I ~/code/microtest/data/atac-seq_SE.bam -O $HOME -S sample3
```

This example is a documented vignette; so just read it and run it to get an idea of how things work.



```python
# %load ../example_pipelines/count_reads.py
#!/usr/bin/env python

"""
Counts reads.
"""

__author__ = "Nathan Sheffield"
__email__ = "nathan@code.databio.org"
__license__ = "GPL3"
__version__ = "0.1"

from argparse import ArgumentParser
import os, re
import sys
import subprocess
import yaml
import pypiper

parser = ArgumentParser(
    description="A pipeline to count the number of reads and file size. Accepts"
    " BAM, fastq, or fastq.gz files.")

# First, add standard arguments from Pypiper.
# groups="pypiper" will add all the arguments that pypiper uses,
# and adding "common" adds arguments for --input and --sample--name
# and "output_parent". You can read more about your options for standard
# arguments in the pypiper docs (section "command-line arguments")
parser = pypiper.add_pypiper_args(parser, groups=["pypiper", "common", "ngs"],
                                    args=["output-parent", "config"],
                                    required=['sample-name', 'output-parent'])

# Add any pipeline-specific arguments if you like here.

args = parser.parse_args()

if not args.input or not args.output_parent:
    parser.print_help()
    raise SystemExit

if args.single_or_paired == "paired":
    args.paired_end = True
else:
    args.paired_end = False

# args for `output_parent` and `sample_name` were added by the standard 
# `add_pypiper_args` function. 
# A good practice is to make an output folder for each sample, housed under
# the parent output folder, like this:
outfolder = os.path.abspath(os.path.join(args.output_parent, args.sample_name))

# Create a PipelineManager object and start the pipeline
pm = pypiper.PipelineManager(name="count",
                             outfolder=outfolder, 
                             args=args)

# NGSTk is a "toolkit" that comes with pypiper, providing some functions
# for dealing with genome sequence data. You can read more about toolkits in the
# documentation

# Create a ngstk object
ngstk = pypiper.NGSTk(pm=pm)

raw_folder = os.path.join(outfolder, "raw/")
fastq_folder = os.path.join(outfolder, "fastq/")

# Merge/Link sample input and Fastq conversion
# These commands merge (if multiple) or link (if single) input files,
# then convert (if necessary, for bam, fastq, or gz format) files to fastq.

# We'll start with a timestamp that will provide a division for this section
# in the log file
pm.timestamp("### Merge/link and fastq conversion: ")

# Now we'll rely on 2 NGSTk functions that can handle inputs of various types
# and convert these to fastq files.

local_input_files = ngstk.merge_or_link(
                        [args.input, args.input2],
                        raw_folder,
                        args.sample_name)

cmd, out_fastq_pre, unaligned_fastq = ngstk.input_to_fastq(
                                            local_input_files,
                                            args.sample_name,
                                            args.paired_end,
                                            fastq_folder)


# Now we'll use another NGSTk function to grab the file size from the input files
#
pm.report_result("File_mb", ngstk.get_file_size(local_input_files))


# And then count the number of reads in the file

n_input_files = len(filter(bool, local_input_files))

raw_reads = sum([int(ngstk.count_reads(input_file, args.paired_end)) 
                for input_file in local_input_files]) / n_input_files

# Finally, we use the report_result() function to print the output and 
# log the key-value pair in the standard stats.tsv file
pm.report_result("Raw_reads", str(raw_reads))

# Cleanup
pm.stop_pipeline()

```
