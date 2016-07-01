
Introduction
=========================

Overview
*************

The target user of Pypiper is a computational scientist comfortable on the command line, who has something like a bash script that would benefit from a layer of "handling code". Pypiper helps you convert that set of shell commands into a production-scale workflow, automatically handling the annoying details (restartablilty, file integrity, logging) to make your pipeline robust and restartable, with minimal learning curve.

Pypiper does not handle any sort of cluster job submission, resource requesting, or parallel dependency management (other than node-threaded parallelism inherent in your commands) -- we use `Looper <http://looper.readthedocs.io/>`_ for that (but you can use whatever you want). Pypiper  just handles a one-sample, sequential pipeline, but now it's robust, restartable, and logged. When coupled with `Looper <http://looper.readthedocs.io/>`_ you get a complete pipeline management system.

Power in simplicity
*********************

Pypiper does not include a new language to learn. You **write your pipeline in python**. Pypiper does not assume you want a complex dependency structure. You write a simple **ordered sequence of commands**, just like a shell script. Pypiper tries to exploit the `Pareto principle <https://en.wikipedia.org/wiki/Pareto_principle>`_ -- you'll get 80% of the features with only 20% of the work of other pipeline management systems.


Installing
*************

The source code lives at Github. You can install directly from GitHub using pip:

.. code-block:: bash

	pip install --user https://github.com/epigen/pypiper/zipball/master


Update with:

.. code-block:: bash

	pip install --user --upgrade https://github.com/epigen/pypiper/zipball/master



Motivation
*************

As I began to put together production-scale pipelines, I found a lot of relevant pipelining systems, but was universally disappointed. For my needs, they were all overly complex. I wanted something simple enough to **quickly write and maintain** a pipeline without having to learn a lot of new functions and conventions, but robust enough to handle requirements like restartability and memory usage monitoring. Everything related was either a pre-packaged pipeline for a defined purpose, or a heavy-duty development environment that was overkill for a simple pipeline. Both of these seemed to be targeted toward ultra-efficient uses, and neither fit my needs: I had a set of commands already in mind -- I just needed a wrapper that could take that code and make it automatically restartable, logged, robust to crashing, easy to debug, and so forth.

If you need a full-blown, datacenter-scale environment that can do everything, look elsewhere. Pypiper's strength is its simplicity. If all you want is a shell-like script, but now with the power of python, and restartability, then Pypiper is for you.

Documentation
*************
Pypiper's documentation is at http://pypiper.readthedocs.org/.

You can also generate docs locally using `sphinx <http://www.sphinx-doc.org/en/stable/install.html>`_. Change your working directory to ``doc`` and run ``make`` to see available documentation formats *e.g.*: ``make html``. The documentation is built in ``doc/build``.

Testing
*************
You can test pypiper by cloning it and running the included unit tests: ``python test_pypiper.py``.

