
# Best practices

Here are some guidelines for how you can design the most effective pipelines.


* **Compartmentalize output into folders**. 
	In your output, keep pipeline steps separate by organizing output into subfolders.

* **Use git for versioning**. 
	If you develop your pipeline in a git repository, Pypiper will automatically record the commit hash when you run a pipeline, making it easy to figure out **exactly** what code version you ran.

* **Record stats as you go**. 
	In other words, don't do all your stats (`report_result()`) and QC at the end; do it along the way. This facilitates monitoring and maximizes availability of statistics even when a pipeline fails.

* **Use looper args**. 
	Even if you're not using looper at first, use `looper_args` and your pipeline will be looper-ready when it comes time to run 500 samples.

* **Use NGSTk early on**. 
	`NGSTk` has lots of useful functions that you will probably need. We've worked hard to make these robust and universal. For example, using NGSTk, you can easily make your pipeline take flexible input formats (FASTQ or BAM). Right now you may always have the same input type (FASTQ, for example), but later you may want your pipeline to be able to work from `bam` files. We've already written simple functions to handle single or multiple BAM or FASTQ inputs; just use this infrastructure (in `NGSTk`) instead of writing your own, and you'll save yourself future headaches.

* **Make some important parameters in the pipeline config, instead of hardcoding them**
	Pypiper makes it painfully easy to use a config file to make your pipeline configurable. Typically you'll start by hard-coding in those parameters in your pipeline steps. But you can select a few important parameters and make them customizable in the pipeline config. Start from the very beginning by making a `yaml` pipeline config file. See an example of a [pipeline config file](configuration.md).
