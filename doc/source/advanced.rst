Advanced
=========================


Toolkits
*************

Pypiper includes optional "toolkits" (right now just one) -- suites of commonly used code snippets which simplifies tasks for a pipeline author. For example, the next-gen sequencing toolkit, NGSTk, which simply provides some convenient helper functions to create common shell commands, like converting from file formats (_e.g._ ``bam_to_fastq()``), merging files (_e.g._ ``merge_bams()``), counting reads, etc. These make it faster to design bioinformatics pipelines in Pypiper, but are entirely optional. Contributions of additional toolkits or functions to an existing toolkit are welcome.


The ``follow`` argument
*************
*Follow functions* let you couple functions to run commands; the functions will only be run if the command is run. This lets you avoid running unnecessary processes repeatedly in the event that you restart your pipeline multiple times (for instance, while debugging later steps in the pipeline).

One of the really useful things in pypiper is that you can pass a python function (we call a "follow function") along in your call to ``run``, and python will call this *follow function* if and only if it runs your command. This is useful for data checks to make sure processes did what you expect. A good example is file conversions; after you convert from one format into another, it's good practice to make sure the conversion completed and you have the same number of lines, or whatever, in the final file. Pypiper's command locking mechanism solves most of this problem, but there are still lots of times you'd like to run a function but only once, right after the command that produced it runs. For example, just counting the number of lines in a file after producing it. You don't need to re-run the counting every time you restart the pipeline, which will skip that step since it is complete. 



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


In this yaml file, the developer of a pipeline records any information the pipeline needs to run that is not related to the Sample being processed. By convension, for consistency across pipelines we use sections called ``tools``, ``resources``, and ``parameters``, but the developer has the freedom to add other sections/variables as needed.

Other information related to a specific run (*e.g.* cpus and memory available) should ideally be passed as command-line arguments.

Pipeline config files by default are named the same as the pipeline with the suffix ``.yaml`` and reside in the same directory as the pipeline code.


Example:

.. code-block:: yaml

	tools:
	  # absolute paths to required tools
	  java:  /home/user/.local/tools /home/user/.local/tools/java
	  trimmomatic:  /home/user/.local/tools/trimmomatic.jar
	  fastqc:  fastqc
	  samtools:  samtools
	  bsmap:  /home/user/.local/tools/bsmap
	  split_reads:  /home/user/.local/tools/split_reads.py  # split_reads.py script; distributed with this pipeline

	resources:
	  # paths to reference genomes, adapter files, and other required shared data
	  resources: /data/groups/lab_bock/shared/resources
	  genomes: /data/groups/lab_bock/shared/resources/genomes/
	  adapters: /data/groups/lab_bock/shared/resources/adapters/

	parameters:
	  # parameters passed to bioinformatic tools, subclassed by tool

	  trimmomatic:
	    quality_encoding: "phred33"
	    threads: 30
	    illuminaclip:
	      adapter_fasta: "/home/user/.local/tools/resources/cpgseq_adapter.fa"
	      seed_mismatches: 2
	      palindrome_clip_threshold: 40
	      simple_clip_threshold: 7
	    slidingwindow:
	      window_size: 4
	      required_quality: 15
	    maxinfo:
	      target_length: 17
	      strictness: 0.5
	    minlen:
	      min_length: 17

	  bsmap:
	    seed_size: 12
	    mismatches_allowed_for_background: 0.10
	    mismatches_allowed_for_left_splitreads: 0.06
	    mismatches_allowed_for_right_splitreads: 0.00
	    equal_best_hits: 100
	    quality_threshold: 15
	    quality_encoding: 33
	    max_number_of_Ns: 3
	    processors: 8
	    random_number_seed: 0
	    map_to_strands: 0






Python process types: Shell vs direct
*************
By default, Pypiper will try to guess what kind of process you want, so for most pipelines, it's probably not necessary to understand the details in this section. However, how you write your commands has some implications for memory tracking, and advanced pipeline authors may want to control the process types that Pypiper uses, so this section covers how these subprocesses work.

Since Pypiper runs all your commands from within python (using the `subprocess` python module), it's nice to be aware of the two types of processes that `subprocess` can handle: **direct processes** and **shell processes**.

**Direct process**: A direct process is one that Python executes directly, from within python. Python retains control over the process completely. For most use cases, you should simply use a direct subprocess (the default) -- this has the advantage of enabling Python to monitor the memory use of the subprocess, because Python retains control over it. This the preferable way of running subprocesses in Python.

**Shell process**: In a shell process, Python first spawns a shell, and then runs the command in that shell. The spawned shell is then directly controlled by Python, but anything the shell does is not; therefore, you lose the ability to monitor memory high water mark because Python does not have direct control over subprocesses run inside a subshell. You must use a shell process if you are using shell operators in your command. For instance, if you use an asterisk (`*`) for wildcard expansion, or a bracket (`>`) for output redirection, or a pipe (`|`) to link processes -- these are commands understood by a shell like Bash, and thus, cannot be run as direct subprocesses in Python.

You can force Pypiper to use one or the other by specifying ``shell=True`` or ``shell=False`` to the ``run`` function. By default Pypiper will try to guess: if your command contains any of the shell process characters ("*", "|", or ">"), it will be run in a shell. Otherwise, it will be run as a direct subprocess.

Harvesting statistics
*************

Pypiper has a neat function called ``get_stat`` that lets you retrieve any value you've reported with ``report_result`` so you could use it to calculate statistics elsewhere in the pipeline. It will retrieve this either from memory, if the calculation of that result happened during the current pipeline run, or from the ``_stats.tsv`` file if the result was reported by an earlier run (or even another pipeline). So you could in theory calculate statistics based on results across pipelines.

An example for how to use this is how we handle calculating the alignment rate in an NGS pipeline:

.. code-block:: python

	x = myngstk.count_mapped_reads(bamfile, args.paired_end)
	pm.report_result("Aligned_reads", x)
	rr = float(pm.get_stat("Raw_reads"))
	pm.report_result("Alignment_rate", round((rr * 100 / float(x), 3))

Here, we use ``get_stat`` to grab a result that we reported previously (with ``report_result``), when we counted the number of ``Raw_reads``. We need this after the alignment to calculate the alignment rate. Later, now that we've reported ``Alignment_rate``, you could harvest this stat again for use with ``pm.get_stat("Alignment_rate")``. This is useful because you could put this block of code in a ``follow`` statement so it may not be executed, but you can still grab a reported result like this even if the execution happened outside of the current pipeline run; you'd only have to do the calculation once.

Any statistics you report like this will be available in summaries made using built-in summary scripts.







