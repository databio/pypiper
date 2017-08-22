
Your first pipeline
***************************

Using pypiper is simple. Your pipeline is a python script, say `pipeline.py`. First, import pypiper, specify an output folder, and create a new ``PipelineManager`` object:

.. code-block:: python

	#!/usr/bin/env python
	import pypiper, os
	outfolder = "pipeline_output/" # Choose a folder for your results
	pm = pypiper.PipelineManager(name="my_pipeline", outfolder=outfolder)

This creates your ``outfolder`` and places a flag called ``my_pipeline_running.flag`` in the folder. It also initializes the log file (``my_pipeline_log.md``) with statistics such as time of starting, compute node, software versions, command-line parameters, etc.

Now, the workhorse of ``PipelineManager`` is the ``run()`` function. Essentially, you just create a shell command as a string in python, and then pass it and its target (a file it creates) to ``run()``. The target is the final output file created by your command. Let's use the built-in ``shuf`` command to create some random numbers and put them in a file called ``outfile.txt``:

.. code-block:: python

	# our command will produce this output file
	target = os.path.join(outfolder, "outfile.txt")
	command = "shuf -i 1-500000000 -n 10000000 > " + target
	pm.run(command, target)


The command (``command``) is the only required argument to ``run()``. You can leave ``target`` empty (pass ``None``). If you **do** specify a target, the command will only be run if the target file does not already exist. If you **do not** specify a target, the command will be run every time the pipeline is run. 

Now string together whatever commands your pipeline requires! At the end, terminate the pipeline so it gets flagged as successfully completed:

.. code-block:: python

	pm.stop_pipeline()

That's it! By running commands through ``run()`` instead of directly in bash, you get a robust, logged, restartable pipeline manager for free!

Go to the next page (:doc:`basic tutorial <tutorials-basic>`) to see a more complicated example.
