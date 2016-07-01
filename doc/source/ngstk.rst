
NGSTk
=========================

An optional feature of pypiper is the accompanying toolkits, such as the next-gen sequencing toolkit, `NGSTk`_, which simply provides some convenient helper functions to create common commands, like converting from file formats (*e.g.* bam to fastq), merging files (*e.g.* merge_bams), counting reads, etc. These make it faster to design bioinformatics pipelines in Pypiper, but are entirely optional.

Example:

.. code-block:: python

	import pypiper
	pm = pypiper.PipelineManager(..., args = args)

	# Create a ngstk object (pass the PipelineManager as an argument)
	ngstk = pypiper.NGSTk(pm = pm)

	# Now you use use ngstk functions
	ngstk.index_bam("sample.bam")


A list of available functions can be found in the :doc:`API <api>` or in the source code for `NGSTk`_.

Contributions of additional toolkits or functions in an existing toolkit are welcome.

.. _NGSTk: https://github.com/epigen/pypiper/blob/master/pypiper/ngstk.py
