
Introduction
=========================

Overview
*************

Pypiper helps you produce pipelines.

Pypiper is a lightweight python toolkit for gluing together restartable, robust
command line pipelines. With Pypiper, simplicity is paramount. Our goal is that
a user could start building useful pipelines using the system in under 15 minutes.
Learning all the possible functions, capabilities, and benefits of Pypiper should
take under an hour. At the same time, using Pypiper should provide immediately clear
and significant advantages over a simple bash script.

The target user of Pypiper is a computational scientist comfortable on the command line, who has something like a bash script that would benefit from a layer of "handling code". Pypiper helps you convert that set of shell commands into a production-scale workflow, automatically handling the annoying details to make your pipeline robust and restartable, with minimal learning curve.

Pypiper does not handle any sort of cluster job submission, resource requesting, or parallel dependency management (other than node-threaded parallelism inherent in your commands). You can use your current setup for those things, and use Pypiper just to produce a robust, restartable, and logged procedural pipeline.

Pypiper provides the following benefits:

.. image:: _static/pypiper.svg

-   **Restartability:** Commands check for their targets and only run if the target needs to be created, much like a `makefile`, making the pipeline pick up where it left off in case it needs to be restarted or extended.
-   **Pipeline integrity protection:** PyPiper uses file locking to ensure that multiple pipeline runs will not interfere with one another -- even if the steps are identical and produce the same files. One run will seemless wait for the other, making it possible to share steps seemlessly across pipelines.
-   **Memory use monitoring:** Processes are polled for high water mark memory use, allowing you to more accurately guage your future memory requirements.
-   **Easy job status monitoring:** Pypiper uses status flag files to make it possible to assess the current state (`running`, `failed`, or `completed`) of hundreds of jobs simultaneously.
-   **Robust error handling:** Pypiper closes pipelines gracefully on interrupt or termination signals, converting the status to `failed`. By default, a process that returns a nonzero value halts the pipeline, unlike in bash, where by default the pipeline would continue using an incomplete or failed result. This behavior can be overridden as desired with a single parameter.
-   **Logging:** Pypiper automatically records the output of your pipeline and its subprocesses, and provides copious information on pipeline initiation, as well as easy timestamping.
-   **Easy result reports:** Pypiper provides functions to put key-value pairs into an easy-to-parse stats file, making it easy to summarize your pipeline results.
-   **Simplicity:** It should only take you 15 minutes to run your first pipeline. The basic documentation is just a few pages long. The codebase itself is also only a few thousand lines of code, making it very lightweight.


Furthermore, Pypiper includes a suite of commonly used pieces of code which the user
use to build it's pipelines.

Installing
*************

.. code-block:: bash

	pip install --user https://github.com/epigen/pypiper/zipball/master


Update with:

.. code-block:: bash

	pip install --user --upgrade https://github.com/epigen/pypiper/zipball/master


Documentation
*************
Pypiper's documentation is at http://pypiper.readthedocs.org/.

You can also generate docs locally using `sphinx <http://www.sphinx-doc.org/en/stable/install.html>`_. With sphinx installed, just change your working directory to `doc` and run `make` to see available documentation formats *e.g.*: `make html`. The documentation will be produced under ``doc/build``.

Toolkits
*************

An optional feature of pypiper is the accompanying toolkits, such as the next-gen sequencing toolkit, NGSTk, which simply provides some convenient helper functions to create common shell commands, like converting from file formats (_e.g._ `bam_to_fastq()`), merging files (_e.g._ `merge_bams()`), counting reads, etc. These make it faster to design bioinformatics pipelines in Pypiper, but are entirely optional. Contributions of additional toolkits or functions to an existing toolkit are welcome.


Motivation
*************

As I began to put together production-scale pipelines, I found a lot of relevant pipelining systems, but was universally disappointed. For my needs, they were all overly complex. I wanted something simple enough to quickly write and maintain a pipeline without having to learn a lot of new functions and conventions, but robust enough to handle requirements like restartability and memory usage monitoring. Everything related was either a pre-packaged pipeline for a defined purpose, or a heavy-duty development environment that was overkill for my needs. Both of these seemed to be targeted toward less experienced developers who sought structure, and neither fit my needs: I had a set of commands already in mind; I just needed a wrapper that could take that code and make it automatically restartable, logged, robust to crashing, easy to debug, and so forth.

If you need a full-blown environment that can do everything, look elsewhere. Pypiper's strength is its simplicity. If all you want is a shell-like script, but now with the power of python, and restartability, then Pypiper is for you.


Testing
*************
You can test pypiper by running ``python test_pypiper.py``, which has some unit tests.

License
*************
Pypiper ___ licensed source code is available at http://github.com/epigen/pypiper/ .

Contributing
*************
We welcome contributions in the form of pull requests; or, if you find a bug or want request a feature, open an issue in https://github.com/epigen/pypiper/issues.

