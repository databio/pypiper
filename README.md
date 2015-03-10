# PyPiper
A lightweight python toolkit for gluing together restartable, robust command line pipelines

# Introduction
---------------
PyPiper is a simple Python module designed to help a computational scientist produce pipelines. As I began to try to put together production-grade pipelines, I found a lot of seeminingly relevant software, but was universally disappointed. I was looking for something simple enough to quickly write and maintain a pipeline, but robust enough to handle requirements like restartability and memory usage monitoring. Everything related was either a pre-packaged pipeline for a defined purpose, or a heavy-duty development environment that was overkill for my needs. Both of these seemed to be targeted toward less experienced developers who sought structure, and neither fit my needs: I had a set of commands already in mind; I just needed a wrapper that could take that code and make it automatically restartable, logged, robust to crashing, easy to debug, and so forth.

The target user of PyPiper is a computational scientist comfortable working on the command line, who has something like a `bash` script that would benefit from a layer of handling code. Converting your shell code into a python script that uses functions provided by PyPiper turns this script into a production-scale workflow that makes it possible to apply the script easily across hundreds or thousands of runs.

PyPiper development is guided by these ideals:
* simplicity - It should only take you 15 minutes to run your first pipeline. The basic documentation is just a few pages long. The codebase itself is also only hundreds of lines of code, making it very lightweight.
* restartability - your pipeline should pick up where it left off if it fails at one step, with minimal manual interaction
* debuggability - It should be easy to identify the source and time of errors.
* transparent - Everything about the run should be logged for furture reference, including pipeline version, compute node, copious timestamps, memory high water mark, etc.

# Your first pipeline
---------------------

The workhorse of PyPiper is the `call_lock()` function. A simple pipeline could look something like this:

```python
import pypiper

# First, generate some random data
cmd = "time shuf -i 1-500000000 -n 100000000 > test.out"
pypiper.call_lock(cmd, "lock.shuf")

cmd = "time paste test.out test.out test.out > pasted.out"
pypiper.call_lock(cmd, "lock.paste")

cmd = "awk 'n < $0 {n=$0}END{print n}' pasted.out"
pypiper.call_lock(cmd, "lock.max")
```




