#!/usr/bin/env python

import pypiper

outfolder = "hello_pypiper_results"  # Choose a folder for your results

# Create a PipelineManager, the workhorse of pypiper
pm = pypiper.PipelineManager(name="hello_pypiper", outfolder=outfolder)

# Timestamps to delineate pipeline sections are easy:
pm.timestamp("Hello!")

# Now build a command-line command however you like, and pass it to pm.run()
target_file = "hello_pypiper_results/output.txt"
cmd = "echo 'Hello, Pypiper!' > " + target_file
pm.run(cmd, target_file)

pm.stop_pipeline()
