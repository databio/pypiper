
Introduction
=========================

Overview
*************

The target user of Pypiper is a computational scientist comfortable on the command line, who has something like a bash script that would benefit from a layer of "handling code". Pypiper helps you convert that set of shell commands into a production-scale workflow, automatically handling the annoying details to make your pipeline robust and restartable, with minimal learning curve.

Pypiper does not handle any sort of cluster job submission, resource requesting, or parallel dependency management (other than node-threaded parallelism inherent in your commands). You can use your current setup for those things, and use Pypiper just to produce a robust, restartable, and logged procedural pipeline.

Installing
*************

The source code lives at Github. You can install directly from GitHub using pip:

.. code-block:: bash

	pip install --user https://github.com/epigen/pypiper/zipball/master


Update with:

.. code-block:: bash

	pip install --user --upgrade https://github.com/epigen/pypiper/zipball/master


Toolkits
*************

Pypiper includes optional "toolkits" -- suites of commonly used code snippets which simplifies tasks for a pipeline author. For example, the next-gen sequencing toolkit, NGSTk, which simply provides some convenient helper functions to create common shell commands, like converting from file formats (_e.g._ ``bam_to_fastq()``), merging files (_e.g._ ``merge_bams()``), counting reads, etc. These make it faster to design bioinformatics pipelines in Pypiper, but are entirely optional. Contributions of additional toolkits or functions to an existing toolkit are welcome.


Motivation
*************

As I began to put together production-scale pipelines, I found a lot of relevant pipelining systems, but was universally disappointed. For my needs, they were all overly complex. I wanted something simple enough to quickly write and maintain a pipeline without having to learn a lot of new functions and conventions, but robust enough to handle requirements like restartability and memory usage monitoring. Everything related was either a pre-packaged pipeline for a defined purpose, or a heavy-duty development environment that was overkill for a simple pipeline. Both of these seemed to be targeted toward less experienced developers who sought structure, and neither fit my needs: I had a set of commands already in mind; I just needed a wrapper that could take that code and make it automatically restartable, logged, robust to crashing, easy to debug, and so forth.

If you need a full-blown, datacenter-scale environment that can do everything, look elsewhere. Pypiper's strength is its simplicity. If all you want is a shell-like script, but now with the power of python, and restartability, then Pypiper is for you.

Documentation
*************
Pypiper's documentation is at http://pypiper.readthedocs.org/.

You can also generate docs locally using `sphinx <http://www.sphinx-doc.org/en/stable/install.html>`_. Change your working directory to ``doc`` and run ``make`` to see available documentation formats *e.g.*: ``make html``. The documentation is built in ``doc/build``.

Testing
*************
You can test pypiper by running ``python test_pypiper.py``, which has some unit tests.

License
*************
Pypiper ___ licensed source code is available at http://github.com/epigen/pypiper/ .

Contributing
*************
We welcome contributions in the form of pull requests; or, if you find a bug or want request a feature, open an issue in https://github.com/epigen/pypiper/issues.

