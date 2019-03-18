# <img src="img/pypiper_logo.svg" class="img-header">

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)

## Introduction


Pypiper is a lightweight python toolkit for gluing together restartable command
line pipelines. With Pypiper, **simplicity is paramount**. It should take less
than 15 minutes to build your first pipeline. Learning all the
[features and benefits](/features) takes just an hour or two. At
the same time, Pypiper provides immediate advantages over a
simple shell script.

Pypiper is an example of a simple  [bioinformatics pipeline framework](
http://databio.org/pipeline_frameworks/). It differs from existing frameworks in its focus on **simplicity** and **sequential pipelines**. 
To employ pypiper, you will just take your bash script and pass those commands through the ``run`` method on a ``PipelineManager`` object. This will give you automatic restartability, process monitoring for memory use and compute time, pipeline status monitoring, copious log output, robust error handling, easy debugging tools, guaranteed file output integrity, and a bunch of useful pipeline development helper functions.




## Installing


Release versions are posted on the GitHub [pypiper releases page](https://github.com/databio/pypiper/releases). You can install the latest release directly from PyPI using pip:

```{console}
pip install --user piper
```

Update `pypiper` with:

```{console}
pip install --user --upgrade piper
```


## Quick start

Build a pipeline `pypiper` via python interface:

```{python}
#!/usr/bin/env python

import pypiper
outfolder = "hello_pypiper_results" # Choose a folder for your results

# Create a PipelineManager, the workhorse of pypiper
pm = pypiper.PipelineManager(name="hello_pypiper", outfolder=outfolder)

# Timestamps to delineate pipeline sections are easy:
pm.timestamp("Hello!")

# Now build a command-line command however you like, and pass it to pm.run()
target_file = "hello_pypiper_results/output.txt"
cmd = "echo 'Hello, Pypiper!' > " + target_file
pm.run(cmd, target_file)

pm.stop_pipeline()
```

Then invoke your pipeline via the command-line:

```{console}
python mypipeline.py --help
```

To begin, check out the [tutorial](tutorial).
