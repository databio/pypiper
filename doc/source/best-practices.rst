
Best practices
=========================

Here are some guidelines for how you can design the most effective pipelines.

* **Parameterize your pipeline** (from the beginning). Pypiper makes it painfully easy to use a config file to make your pipeline configurable. Start from the very beginning by making a ``yaml`` pipeline config file.

* **Use flexible inputs**. Right now you may always have the same input type (fastq, for example), but later you may want your pipeline to be able to work from ``bam`` files. We've already written simple functions to handle single or multiple bam or fastq inputs; just use this infrastructure (in NGSTk) instead of writing your own, and you'll save yourself future headaches.

* **Use looper args**. Even if you're not using looper at first, use ``looper_args`` and your pipeline will be looper-ready when it comes time to run 500 samples.

* **Compartmentalize with folders**. In your output, keep pipeline steps separate by organizing output into folders.

* **Use git for versioning**. If you develop your pipeline in a git repository, Pypiper will automatically record what commit you run, making it easy to figure out what code version you ran.

* **Record stats as you go**. In other words, don't do all your stats and QC at the end, do it along the way. This makes it easy for you to monitor pipeline performance, and couples stats with how far the pipeline makes it, so you could make use of a partially completed (or even ultimately failed) pipelines.
