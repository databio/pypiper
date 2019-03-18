
# What are toolkits?

Pypiper provides a limited set of functions that are all quite generic; they simple accept command-line commands and run them. You could use this to produce a pipeline in any domain.

To add to this, you may be interested in building some of your own convenience functions that make it easier for you to piece together commands. It's really easy to  create your own library of python functions by creating a python package. Then, you just need to import your package in your pipeline script and make use of the common functions.

We refer to this type of package as a "toolkit", and Pypiper also includes an built-in toolkit called NGSTk (next-generation sequencing toolkit). NGSTk simply provides some convenient helper functions to create common shell commands, like converting from file formats (_e.g._ `bam_to_fastq()`), merging files (_e.g._ `merge_bams()`), counting reads, etc. These make it faster to design bioinformatics pipelines in Pypiper, but are entirely optional. Contributions of additional toolkits or functions to an existing toolkit are welcome.

More details about the NGSTk toolkit can be found on the next page under [NGSTk](ngstk).
