
Using pipelines
=========================

Pypiper pipelines are python scripts. There is no special requirement or syntax, you run it as you would any python script:

.. code-block:: bash

	python pipeline.py

Or, if you make the script executable (``chmod o+x pipeline.py``) and it has a shebang at the top (``#!/usr/bin/env python``), you can execute it directly:

.. code-block:: bash

	./pipeline.py


Now, you'll need to figure out the command-line arguments required by the pipeline. Usage will vary based on the script, and the pipeline author determines what command line arguments their pipeline will recognize. You should look in the documentation for your pipeline, or you can often figure out the command-line interface by passing the ``--help`` argument, like so: 

.. code-block:: bash

	pipeline.py --help


Universal pypiper options
^^^^^^^^^^^^^^^^^^^^^^^^^^

With that said, there are a few universal (Pypiper-added) options that are frequently (but not necessarily always) honored by pypiper pipelines. These default pypiper arguments are detailed below:

  - ``-R, --recover``
  	Recover mode, overwrite locks. This argument will tell pypiper to recover from a failed previous run. Pypiper will execute commands until it encounters a locked file, at which point it will re-execute the failed command and continue from there.

  - ``-F, --follow``
  	Force run follow-functions. By default, follow-functions are only run if their corresponding ``run`` command was run; with this option you can force all follow functions to run. This is useful for regenerating QC data on existing output. For more details, see :ref:`the follow argument <the_follow_argument>`.

  - ``-D, --dirty``
  	Make all cleanups manual. By default, pypiper pipelines will delete any intermediate files. For debugging, you may want to turn this option off -- you can do that by specifying **dirty mode**.

  - ``-N, --new-start``
  	New start mode. This flag will tell pypiper to start over, and run every command, even if its target output already exists.
