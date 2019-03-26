# <img src="img/pypiper_logo.svg" class="img-header">

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)

## Introduction

Pypiper is a **development-oriented** pipeline framework. Pypiper pipelines are:

1. written in pure python, so they require learning no new language;
2. simple to update and maintain, so they respond well to changes;
3. simple to understand for an outsider, so they can be approached by others.


These traits make pypiper ideally suited for **pipelines under active development**.

With Pypiper, **simplicity is paramount**. Prerequisites are few: base python and 2 common packages (`pyyaml` and `psutil`). It should take fewer than 15 minutes to build your first pipeline and only an hour or two to learn the advanced features.
Pypiper provides automatic restartability, process monitoring for time and memory use, status monitoring, copious log output, robust error handling, easy debugging tools, guaranteed file output integrity, and a bunch of useful pipeline development helper functions.  Read more about the [pypiper philosophy](philosophy).

## Installing

Release versions are posted on the GitHub [pypiper releases page](https://github.com/databio/pypiper/releases). You can install the latest release directly from PyPI using `pip`.

Global scope for single user:
```{console}
pip install --user --upgrade piper
```

Within an active virtual environment:
```{console}
pip install --upgrade piper
```

## Quick start

To employ pypiper, you build something like a shell script, but pass the commands through the `run` method on a `PipelineManager` object. Build your pipeline in **pure python**:

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

Pypiper differs from existing frameworks in its focus on **simplicity**. Pypiper requires learning no new language, as **pipelines are written in pure python**. Pypiper is geared toward **developing pipelines** that are contained in a single file, easy to update, and easy to understand. Read more about
