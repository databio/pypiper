
Your first pipeline
***************************

Using pypiper is simple. First, import pypiper, specify an output folder, and create a new ``PipelineManager`` object:

.. code-block:: python

	import pypiper, os
	outfolder = "pipeline_output/" # Choose a folder for your results
	pipeline = pypiper.PipelineManager(name="my_pipeline", outfolder=outfolder)

This creates your ``outfolder`` and places a flag called ``my_pipeline_running.flag`` in the folder. It also initializes the log file (``my_pipeline_log.md``) with statistics such as time of starting, compute node, software versions, command-line parameters, etc.

Now, the workhorse of ``PipelineManager`` is the ``run()`` function. Essentially, you just create a shell command as a string in python, and then pass it and its target (a file it creates) to ``run()``. The target is the final output file created by your command.

.. code-block:: python

	# our command will produce this output file
	target = os.path.join(outfolder, "outfile.txt")
	cmd = "shuf -i 1-500000000 -n 10000000 > " + target
	pipeline.run(command, target)


The command (``cmd``) is the only required argument to ``run()``. You can leave ``target`` empty (pass ``None``). If you **do** specify a target, the command will only be run if the target file does not already exist. If you **do not** specify a target, the command will be run every time the pipeline is run. 

Now string together whatever commands your pipeline requires! At the end, terminate the pipeline so it gets flagged as successfully completed:

.. code-block:: python

	pipeline.stop_pipeline()

That's it! By running commands through ``run()`` instead of directly in bash, you get a robust, logged, restartable pipeline manager for free!

To see an example of a simple pipline, look in the `example_pipelines` folder in this respository (also listed here under tutorials), which are thoroughly commented to act as vignettes. This is the best way to learn how to use Pyipiper.
