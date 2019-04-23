# FAQ

## How can I run my pipeline on more than 1 sample?

Pypiper only handles individual-sample pipelines. To run it on multiple samples, write a loop, or use [looper](http://looper.readthedocs.io/). Dividing multi-sample handling from individual sample handling is a conceptual advantage that allows us to write a nice, universal, generic sample-handler that you only have to learn once.

## What cluster resources can pypiper use?

Pypiper is compute-agnostic. You run it wherever you want. If you want a nice way to submit pipelines for samples any cluster manager, check out [looper](http://looper.readthedocs.io/), which can run your pipeline on any compute infrastructure using the [divvy python package](http://code.databio.org/divvy).

## What does it mean for a sample to be in the "waiting" state?

Waiting means `pypiper` encountered a file lock, but no recovery flag. So the pipeline thinks a process (from another run or another process) is currently writing that file. It periodically checks for the lock file to disappear, and assumes that the other process will unlock the file when finished. If you are sure there's not another process writing to that file, you can get `pypiper` to continue by deleting the corresponding `lock` file. In the future, you can use `pypiper's` recover mode (`-R`) to automatically restart a process when a `lock` file is found, instead of waiting.

## What is the 'elapsed time' in output?

The "elapsed" time is referring to the amount of time since the preceding timestamp, not since the start of the pipeline. Timestamps are all displayed with a flag: `_TIME_`. The total cumulative time for the pipeline is displayed only at the end.

## How should I run a QC step to check results of one of my commands?

Usually, you only want to run a QC step if the result was created in the same pipeline run. There's no need to re-run that step if you have to restart the pipeline due to an error later on. If you use `run()` for these steps, then they'll need to run each time the pipeline runs. Instead, this is exactly why we created [the follow argument](../advanced-run-method/#the-follow-argument) This option lets you couple a QC step to a `run()` call, so it only gets excecuted when it is required.

## How do I solve installation errors involving `psutil` and/or a compiler like `gcc` or `clang`?

If you have trouble with installation and it looks like one of these pieces of software is involved, please check the [`psutil` installation guide](https://github.com/giampaolo/psutil/blob/master/INSTALL.rst).

