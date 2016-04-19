Advanced
=========================

Pipeline command-line arguments
*************
To take full advantage of Pypiper (make your pipeline recoverable, etc.), you need to add command-line options to the ``PipelineManager``. Use the typical Python `argparse module <https://docs.python.org/2/library/argparse.html>`_,  add Pypiper args to it (with ``add_pypiper_args()``), and then pass this to your `PipelineManager`.

.. code-block:: python

	import pypiper, os, argparse
	parser = ArgumentParser(description='Pipeline')

	# add any custom args here
	parser = pypiper.add_pypiper_args(parser)
	

	outfolder = "pipeline_output/" # Choose a folder for your results

	# include your parser as you construct a PipelineManager
	pipeline = pypiper.PipelineManager(name="my_pipeline", outfolder=outfolder, \
						args = parser)

Your pipeline will then enable a few **pypiper keyword arguments**, which are: ``recover``, ``fresh``, and ``dirty``. As a side bonus, all arguments (including any of your custom arguments) will be recorded in the log outputs. What do these keyword arguments do?

- **recover**: for a failed pipeline run, start off at the last succesful step. Overwrites the ``failed`` flag.
- **fresh**: NOT IMPLEMENTED! The plan is for this to just recreate everything, even if it exists.
- **dirty**: Disables automatic cleaning of temporary files, so all intermediate files will still exist after a pipeline run (either sucessful or failed). Useful for debugging a pipeline even if it succeeds.

Additional stuff (to be explained more thorougly soon): You can also add some other advanced parameters by tweaking the arguments pased to ``add_pypiper_args()``:.

 ``config_file``, ``output_parent``, ``cores``, and ``mem``. You would then be able to access these parameters from the pipeline manager, for example using ``pipeline.cores`` or ``pipeline.mem`` in your pipeline.

The most significant of these special keywords is the ``config_file`` argument, which leads us to the concept of ``pipeline config files``:

Pipeline config files
*************
Optionally, you may choose to conform to our standard for parameterizing pipelines, which enables you to use some powerful features of the pypiper system.

If you write a pipeline config file in ``yaml`` format and name it the same thing as the pipeline (but ending in ``.yaml`` instead of ``.py``), pypiper will automatically load and provide access to these configuration options, and make it possible to pass customized config files on the command line. This is very useful for tweaking a pipeline for a similar project with slightly different parameters, without having to re-write the pipeline.

It's easy: just load the ``PipelineManager`` with ``args`` (as above), and you have access to the config file automatically in ``pipeline.config``.

For example, in ``myscript.py`` you write:

.. code-block:: python

	pipeline = pypiper.PipelineManager(name="my_pipeline", outfolder=outfolder, \
						args = parser)


And in the same folder, you include a ``yaml`` called ``myscript.yaml``:

.. code-block:: yaml

	settings:
	  setting1: True
	  setting2: 15

Then you can access these settings automatically in your script using:

.. code-block:: python

	pipeline.config.settings.setting1
	pipeline.config.settings.setting2


Python process types: Shell vs direct
*************
Since Pypiper runs all your commands from within python (using the `subprocess` python module), you need to be aware of the two types of processes that `subprocess` can handle: Shell process and direct processes. Choosing between the two is as simple as specifying `shell=True` (the default is a direct subprocess) to `run()`, which is directly passed on as the `shell` parameter to `subprocess.Popen()`. The difference is that Python runs the process either as a subprocess directly from within python, or it first spawns a shell, and then runs the subprocess in the shell. When should you use each?

**Direct process**: For most use cases, you should simply use a direct subprocess (the default) -- this has the advantage of enabling Python to monitor the memory use of the subprocess, because Python retains control over it. This is considered by the community to be the preferable way of running subprocesses in Python.

**Shell process**: You must use a `shell=True` process if you are using shell operators in your command. For instance, if you use an asterisk (`*`) for wildcard expansion, or a bracket (`>`) for output redirection, or a pipe (`|`) to link processes -- these are commands understood by a shell like Bash, and thus, cannot be run as direct subprocesses by Python, but must be run in a shell. `Subprocess` (and therefore `Pypiper`) can understand these commands just fine, but you lose the ability to monitor memory high water mark because Python does not have direct control over subprocesses run inside a subshell.

