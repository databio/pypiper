# PyPiper
A lightweight python toolkit for gluing together restartable, robust command line pipelines

Introduction
---------------
PyPiper helps you produce pipelines.

The target user of PyPiper is a computational scientist comfortable working on the command line, who has something like a `bash` script that would benefit from a layer of handling code. PyPiper will enable you to quickly convert that set of raw commands into a production-scale workflow, automatically handling the annoying details to make your pipeline robust and deployable.

Benefits of using PyPiper
-------------------------
* Restartability - Commands check for their targets and only run if the target needs to be created, much like a `makefile`, making the pipeline pick up where it left off in case it needs to be restarted or extended.
* Pipeline integrity protection - PyPiper uses file locking to ensure that multiple pipeline runs will not interfere with one another -- even if the steps are identical and produce the same files. One run will seemless wait for the other, making it possible to share steps seemlessly across pipelines.
* Memory use monitoring - Processes are polled for high water mark memory use, allowing you to more accurately guage your future memory requirements.
* Easy job status monitoring - PyPiper uses status flags to make it possible to assess the current state (`running`, `failed`, or `completed`) of hundreds of jobs simultaneously.
* Robust error handling - PyPiper closes pipelines gracefully on interrupt or termination signals, converting the status to `failed`.
* Logging - PyPiper automatically records the output of your pipeline and its subprocesses, and provides copious information on pipeline initiation, as well as easy timestamping.
* Easy result reports - PyPiper provides functions to put key-value pairs into an easy-to-parse stats file, making it easy to summarize your pipeline results.
* Simplicity - It should only take you 15 minutes to run your first pipeline. The basic documentation is just a few pages long. The codebase itself is also only hundreds of lines of code, making it very lightweight.


# Your first pipeline
---------------------

The workhorse of PyPiper is the `call_lock()` function. To see an example of a simple pipline, look at the `sample_pipeline.py` script in this repository.




Motivation
----------
As I began to try to put together production-grade pipelines, I found a lot of seeminingly relevant software, but was universally disappointed. I was looking for something simple enough to quickly write and maintain a pipeline, but robust enough to handle requirements like restartability and memory usage monitoring. Everything related was either a pre-packaged pipeline for a defined purpose, or a heavy-duty development environment that was overkill for my needs. Both of these seemed to be targeted toward less experienced developers who sought structure, and neither fit my needs: I had a set of commands already in mind; I just needed a wrapper that could take that code and make it automatically restartable, logged, robust to crashing, easy to debug, and so forth.


