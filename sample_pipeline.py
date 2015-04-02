#!/usr/bin/python2.7
"""Getting Started: A simple sample pipeline built using pypiper."""

# First, make sure you can import the pypiper package in some way
# Here I just add the path to a cloned repo before importing.

import os

pypiper_dir = "/home/nsheffield/uncorrupted/pypiper"
os.sys.path.insert(0,pypiper_dir)

import pypiper

# Import any supplemental modules, if you'd like
from pypiper import ngstk


# Create a Pypiper instance (don't forget to name it!),
# and then start the pipeline!

mypiper = pypiper.Pypiper(name="sample_pipeline", outfolder="pipeline_output/")
mypiper.start_pipeline()

# Now just build command strings, and use the call_lock function
# to execute them in order.


# First, generate some random data
# You must use shell=True here because of redirection (>), which
# is a shell process, and can't be run as a python subprocess.
cmd = "shuf -i 1-500000000 -n 100000000 > pipeline_output/test.out"
mypiper.call_lock(cmd, target="pipeline_output/test.out", shell=True)

# Now copy the data into a new file.
# No pipes or redirects, so this does not require shell function.
# (this should be the most common use case).
cmd = "cp pipeline_output/test.out pipeline_output/copied.out"
mypiper.call_lock(cmd, target="pipeline_output/copied.out")

# You can also string multiple commands together, which will execute
# in order as a group to create the final target.
cmd1 = "sleep 5"
cmd2 = "touch pipeline_output/touched.out"
mypiper.call_lock([cmd1, cmd2], target="pipeline_output/touched.out")

# A command without a target will run every time.
# Find the biggest line
#cmd = "awk 'n < $0 {n=$0} END{print n}' pipeline_output/test.out"
#mypiper.call_lock(cmd, "lock.max", shell=True)

# Use report_result to print and log key-value pairs in the stats file:
import subprocess
last_entry = subprocess.check_output("tail -n 1 pipeline_output/copied.out", shell=True)
mypiper.report_result("last_entry", last_entry)


# Now, stop the pipeline to complete gracefully.
mypiper.stop_pipeline()

