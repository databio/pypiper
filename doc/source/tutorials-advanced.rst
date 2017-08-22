Advanced tutorial
********************************************************************************

Here we have a more advanced bioinformatics pipeline that adds some new concepts. This is a simple script that takes an input file and returns the file size and the number of sequencing reads in that file. This example uses a function from from the built-in :doc:`NGSTk toolkit <ngstk>`. In particular, this toolkit contains a few handy functions that make it easy for a pipeline to accept inputs of various types. So, this pipeline can count the number of reads from files in ``BAM`` format, or ``fastq`` format, or ``fastq.gz`` format. You can also use the same functions from NGSTk to develop a pipeline to do more complicated things, and handle input of any of these types.

First, grab this pipeline. Download `count_reads.py <https://github.com/epigen/pypiper/blob/master/example_pipelines/count_reads.py>`_, make it executable (``chmod 755 count_reads.py``), and then run it with ``./count_reads.py``). 

You can grab a few small data files in the `microtest repository <https://github.com/epigen/microtest/tree/master/config>`_. Run a few of these files like this:


.. code-block:: shell

	./count_reads.py -I ~/code/microtest/data/rrbs_PE_R1.fastq.gz -O $HOME -S sample1
	./count_reads.py -I ~/code/microtest/data/rrbs_PE_fq_R1.fastq -O $HOME -S sample2
	./count_reads.py -I ~/code/microtest/data/atac-seq_SE.bam -O $HOME -S sample3

This example is a documented vignette; so just read it and run it to get an idea of how things work.

.. literalinclude:: ../../example_pipelines/count_reads.py


