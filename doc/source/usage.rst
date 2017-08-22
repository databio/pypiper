
Using pipelines
=========================

Pypiper pipelines are python scripts. There is no special requirement or syntax. To run a pypiper pipeline, you run it as you would any python script. Usage will vary based on the script, and the pipeline author determines what command line arguments their pipeline will recognize.

You can often figure out the command-line interface by passing the ``--help argument``, like so: ``script.py --help``. 

Many pipelines employ default pypiper arguments which are detailed below:

  - ``-R, --recover``
  	Recover mode, overwrite locks. This argument will tell pypiper to recover from a failed previous run. Pypiper will execute commands until it encounters a locked file, at which point it will re-execute the failed command and continue from there.

  - ``-F, --follow``
  	Force run follow-functions. For more details, see :ref:`the follow argument <the_follow_argument>`.
