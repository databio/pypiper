
Your first pipeline
***************************

Using pypiper is simple. First, import pypiper, specify an output folder, and create a new ``PipelineManager`` object:

.. code-block:: python

	import pypiper, os
	outfolder = "pipeline_output/" # Choose a folder for your results
	pipeline = pypiper.PipelineManager(name="my_pipeline", outfolder=outfolder)

Now, the workhorse of ``PipelineManager`` is the ``run()`` function. Essentially, you just create a shell command as a string in python, and then pass it and its target (a file it creates) to ``run()``. 

.. code-block:: python

	target = os.path.join(outfolder, "outfile.txt")
	command cmd = "shuf -i 1-500000000 -n 10000000 > " + target
	pipeline.run(command, target, shell=True)


The target is the final output file created by your command. You can also leave it empty (pass ``None``) if there is no output. Now string together whatever commands your pipeline requires!

At the end, terminate the pipeline so it gets flagged as successfully completed:

.. code-block:: python

	pipeline.stop_pipeline()

That's it! By running commands through ``run()`` instead of directly in bash, you get a robust, logged, restartable pipeline manager for free!

There's more information about ``shell=True`` later.

To see an example of a simple pipline, look in the `example_pipelines` folder in this respository, which are thoroughly commented to act as vignettes. This is the best way to learn how to use Pyipiper.

