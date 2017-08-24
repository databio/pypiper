Pipeline config files
================================================================================

.. _pipeline_config_files:


Optionally, you may choose to conform to our standard for parameterizing pipelines, which enables you to use some powerful features of the pypiper system.

If you write a pipeline config file in ``yaml`` format and name it the same thing as the pipeline (but ending in ``.yaml`` instead of ``.py``), pypiper will automatically load and provide access to these configuration options, and make it possible to pass customized config files on the command line. This is very useful for tweaking a pipeline for a similar project with slightly different parameters, without having to re-write the pipeline.

It's easy: just load the ``PipelineManager`` with ``args`` (as described in :doc:`command-line arguments <develop-arguments>`), and you have access to the config file automatically in ``pipeline.config``.

For example, in ``myscript.py`` you write:

.. code-block:: python

	parser = pypiper.add_pipeline_args(parser, args=["config"])
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


This ``yaml`` file is useful for any settings the pipeline needs that is not related to the input Sample (which should be passed on the command-line). By convention, for consistency across pipelines, we use sections called ``tools``, ``resources``, and ``parameters``, but the developer has the freedom to add other sections/variables as needed.

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



