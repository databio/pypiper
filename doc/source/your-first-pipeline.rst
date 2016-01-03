
Your first pipeline
=========================

Using Pypiper is simple. First, import pypiper, and then create a new PipelineManager object:

.. code-block:: python

	import pypiper, os
	outfolder = "pipeline_output/" # Choose a folder for your results
	pipeline = pypiper.PipelineManager(name="my_pipeline", outfolder=outfolder)

Now, the workhorse of PipelineManager is the `run()` function. Essentially, you just create a shell command as a string in python, and then pass it and its output to ``run()``. 

.. code-block:: python

	target = os.path.join(outfolder, "outfile.txt")
	command cmd = "shuf -i 1-500000000 -n 10000000 > " + target
	pipeline.run(command, target, shell=True)


The target is the final output file created by your command. You can also leave it empty (pass ``None``) if there is no output. Now string together whatever commands your pipeline requires, and then at the end, terminate the pipeline so it gets flagged as successfully completed:

.. code-block:: python

	pipeline.stop_pipeline()

That's it! By running this command through ``run()`` instead of directly in bash, you get a robust, logged, restartable pipeline manager for free! There's more information about ``shell=True`` processes below.
