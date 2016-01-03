
NGSTk
=========================

An optional feature of pypiper is the accompanying toolkits, such as the next-gen sequencing toolkit, `NGSTk`_, which simply provides some convenient helper functions to create common commands, like converting from file formats (*e.g.* bam to fastq), merging files (*e.g.* merge_bams), counting reads, etc. These make it faster to design bioinformatics pipelines in Pypiper, but are entirely optional.

Example:

.. code-block:: python

	from pypiper.ngstk import NGSTk

	tk = NGSTk()
	tk.index_bam("sample.bam")

Contributions of additional toolkits or functions in an existing toolkit are welcome.

.. _NGSTk: https://github.com/epigen/pypiper/pypiper/ngstk.py
