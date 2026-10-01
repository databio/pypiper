
# NGSTk - Next Gen Sequencing Toolkit

Pypiper functions are generic; they simply accept command-line commands and run them. You could use this to produce a pipeline in any domain. To add to this, it's helpful to build convenience functions specific to your scientific domain. It's really easy to create your own library of python functions by creating a python package. Then, you just need to import your package in your pipeline script and make use of the common functions. We refer to this type of package as a "toolkit".

Pypiper includes a built-in toolkit called NGSTk (next-generation sequencing toolkit). NGSTk simply provides some convenient helper functions to create common shell commands, like converting from file formats (_e.g._ `bam_to_fastq()`), merging files (_e.g._ `merge_bams()`), counting reads, etc. These make it faster to design bioinformatics pipelines in Pypiper, but are entirely optional.

Here's how to use `NGSTk`:

```{python}
import pypiper
pm = pypiper.PipelineManager(..., args = args)

# Create a ngstk object (pass the PipelineManager as an argument)
ngstk = pypiper.NGSTk(pm = pm)

# Now you use use ngstk functions
cmd = ngstk.index_bam("sample.bam")
pm.run(cmd, target="sample.bam")
```

A complete list of functions is in the [API](code/ngstk.md) or in the [source code for NGSTk](https://github.com/databio/pypiper/blob/master/pypiper/ngstk.py).
