
Using Pipelines
=========================

Pypiper pipelines are python scripts. There is no special requirement or syntax. To run a pypiper pipeline, you run it as you would any basic script, usage will vary based on the script. Pipeline authors determine what command line arguments their pipeline will recognize. You can usually figure this out by passing the ``--help argument``, like so: ``script.py --help``. You should check the individual Pypiper script you're using to see what it does.

Many pipelines employ default pypiper arguments which are detailed below:

  - ``-R, --recover``
  	Recover mode, overwrite locks. This argument will tell pypiper to recover from a failed previous run. Pypiper will execute commands until it encounters a locked file, at which point it will re-execute the failed command and continue from there.
  - ``-F, --follow``
  	Run follow-functions. :ref:`the follow argument <The-follow-argument>`.
