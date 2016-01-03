
Pypiper documentation
=========================

Pypiper is a lightweight python toolkit for gluing together restartable, robust
command line pipelines. With Pypiper, simplicity is paramount. Our goal is that
a user could start building useful pipelines using the system in under 15 minutes.
Learning all the possible functions, capabilities, and benefits of Pypiper should
take under an hour. At the same time, using Pypiper should provide immediately clear
and significant advantages over a simple bash script.

Pypiper provides the following virtues:

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

Documentation
*************
Pypiper's documentation is at http://pypiper.readthedocs.org/.

You can generate the documentation yourself, using ``make``. Just change your working directory to `doc` and run `make` to see available documentation formats *e.g.*: `make html`. The documentation will be produced under ``doc/build``.

Testing
*************
You can test pypiper by running ``python test_pypiper.py``, which has some unit tests.

Licence
*************
Pypiper ___ licensed source code is available at http://github.com/epigen/pypiper/ .

Contributing
*************
We welcome contributions in the form of pull requests or if you find a bug or want request a feature, open an issue in https://github.com/epigen/pypiper/issues .


Contents
^^^^^^^^

.. toctree::
   :maxdepth: 2

   your-first-pipeline.rst
   tutorials.rst
   ngstk.rst
   api.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
