# <img src="img/pypiper_logo.svg" class="img-header">

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)

## Introduction

Pypiper is a **development-oriented** pipeline development framework. With Pypiper, **simplicity is paramount**. It should take fewer than 15 minutes to build your first pipeline and only an hour or two to learn the advanced features. Pypiper pipelines **are simple to update and maintain**, making it **geared toward pipelines under active developement**.

Pypiper is an example of a simple [bioinformatics pipeline framework](
http://databio.org/pipeline_frameworks/). 
To employ pypiper, you will just take your bash script and pass those commands through the ``run`` method on a ``PipelineManager`` object. This will give you automatic restartability, process monitoring for memory use and compute time, pipeline status monitoring, copious log output, robust error handling, easy debugging tools, guaranteed file output integrity, and a bunch of useful pipeline development helper functions.

## Installing

Release versions are posted on the GitHub [pypiper releases page](https://github.com/databio/pypiper/releases). You can install the latest release directly from PyPI using pip:

```{console}
pip install --user --upgrade piper
```

## Quick start

Build your pipeline in **pure python**:

```{python}
#!/usr/bin/env python

import pypiper
outfolder = "hello_pypiper_results" # Choose a folder for your results

# Create a PipelineManager, the workhorse of pypiper
pm = pypiper.PipelineManager(name="hello_pypiper", outfolder=outfolder)

# Timestamps to delineate pipeline sections are easy:
pm.timestamp("Hello!")

# Now build a command and pass it to pm.run()
target_file = "hello_pypiper_results/output.txt"
command = "echo 'Hello, Pypiper!' > " + target_file
pm.run(command, target_file)

pm.stop_pipeline()
```

Then invoke your pipeline via the command-line:

```{console}
python my_pipeline.py --help
```

## Pypiper strengths

Pypiper differs from existing frameworks in its focus on **simplicity**. Pypiper requires learning no new language, as **pipelines are written in pure python**. Pypiper is geared toward **sequential pipelines** that are contained in a single file, easy to update, and easy to understand. Read more about the [pypiper philosophy](philosophy).