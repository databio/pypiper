#!/usr/bin/env python
import pypiper
outfolder = "hello_pypiper_results" # Choose a folder for your results
pm = pypiper.PipelineManager(name="hello_pypiper", outfolder=outfolder)


pm.timestamp("Hello!")
target_file = "hello_pypiper_results/output.txt"
cmd = "echo 'Hello, Pypiper!' > " + target_file
pm.run(cmd, target_file)

pm.stop_pipeline()