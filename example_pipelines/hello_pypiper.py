#!/usr/bin/env python

import pypiper

with pypiper.PipelineManager("hello_pypiper", "hello_pypiper_results") as pm:
    pm.timestamp("Hello!")
    target_file = "hello_pypiper_results/output.txt"
    cmd = f"echo 'Hello, Pypiper!' > {target_file}"
    pm.run(cmd, target_file)
