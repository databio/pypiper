# <img src="img/pypiper_logo.svg" class="img-header">

<p align="center">
<a href="https://pep.databio.org" alt="PEP compatible"><img src="https://pepkit.github.io/img/PEP-compatible-green.svg"/></a>
<a href="https://github.com/databio/pypiper/actions/workflows/run-pytest.yml" alt="Run pytests"><img src="https://github.com/databio/pypiper/workflows/Run%20pytests/badge.svg"/></a>
<a href="https://pypi.org/project/piper" alt="PyPI badge"><img src="https://img.shields.io/pypi/v/piper"/></a>
<a href="https://github.com/databio/pypiper" alt="GitHub source code"><img src="https://img.shields.io/badge/source-github-354a75?logo=github"/></a>
</p>



## What is pypiper?

`Pypiper` is a **development-oriented** pipeline framework. It is a Python package that helps you write robust pipelines directly in Python, handling mundane tasks like restartability, monitoring for time and memory use, monitoring job status, copious log output, robust error handling, easy debugging tools, and guaranteed file output integrity.



## What makes pypiper better? 
With Pypiper, **simplicity is paramount**. Prerequisites are few: base Python and 2 common packages (`pyyaml` and `psutil`). It should take fewer than 15 minutes to build your first pipeline and only an hour or two to learn the advanced features. Pypiper pipelines are:

1. written in pure Python, so they do not require learning a new language;
2. easy to modify, so they are simple to update and maintain;
3. simple to understand for an outsider, so they can be approached by others.

These traits make pypiper ideally suited for **pipelines under active development**. Read more about the [pypiper philosophy](philosophy.md).

## Installing

Releases are posted as [GitHub releases](https://github.com/databio/pypiper/releases), or you can install from PyPI using `pip`:

Global scope for single user:
```console
pip install --user --upgrade piper
```

Within an active virtual environment:
```console
pip install --upgrade piper
```

## Quick start

To employ pypiper, you build something like a shell script, but pass the commands through the `run` method on a `PipelineManager` object. Build your pipeline in **pure python**:

```python
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

```console
python my_pipeline.py --help
```

## Pypiper strengths

Pypiper differs from existing frameworks in its focus on **simplicity**. Pypiper requires learning no new language, as **pipelines are written in pure python**. Pypiper is geared toward **developing pipelines** that are contained in a single file, easy to update, and easy to understand.
