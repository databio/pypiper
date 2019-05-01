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
