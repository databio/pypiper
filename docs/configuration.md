# Pipeline configuration files

If you write a pipeline config file in `yaml` format and name it the same thing as the pipeline (but replacing `.py` with `.yaml`), pypiper will automatically load and provide access to these configuration options, and make it possible to pass customized config files on the command line. This is very useful for tweaking a pipeline for a similar project with slightly different parameters, without having to re-write the pipeline.

It's easy: just load the `PipelineManager` with `args` (as described in [command-line arguments](cli.md)), and you have access to the config file automatically in in `pipeline.config`.

For example, in `myscript.py` you write:

```{python}
parser = pypiper.add_pipeline_args(parser, args=["config"])
pipeline = pypiper.PipelineManager(name="my_pipeline", outfolder=outfolder, \
					args = parser)
```

And in the same folder, you include `myscript.yaml`:



	my_section:
	  setting1: True
	  setting2: 15

Then you can access these settings automatically in your script using:



	pipeline.config.my_section.setting1
	pipeline.config.my_section.setting2


This `yaml` file is useful for any parameters *not related to the input Sample* (which should be passed on the command-line). By convention, for consistency across pipelines, we use sections called `tools`, `resources`, and `parameters`, but the developer has the freedom to add other sections/variables as needed.

Here's a more realist example pipeline configuration file:


```{yaml}
# paths to required tools
tools:
  java:  "/home/user/.local/tools/java"
  trimmomatic:  "/home/user/.local/tools/trimmomatic.jar"
  fastqc:  "fastqc"
  samtools:  "samtools"
  bsmap:  "/home/user/.local/tools/bsmap"
  split_reads:  "/home/user/.local/tools/split_reads.py"

# paths to reference genomes, adapter files, and other required shared data
resources:
  resources: "/data/groups/lab_bock/shared/resources"
  genomes: "/data/groups/lab_bock/shared/resources/genomes/"
  adapters: "/data/groups/lab_bock/shared/resources/adapters/"

# parameters passed to bioinformatic tools, subclassed by tool
parameters:
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
```


