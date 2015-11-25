# Command line args

You may choose to use Pypiper to add arguments (via `ArgumentParser`) to your pipeline. There are several sets of options, which may be mixed and matched. You may add all of the options with `add_pypiper_args(all_args=True)`


## Pypiper arguments
These are the default options which will be added by `add_pypiper_args()`. They are used directly by the PipelineManager object and should probably be included with every pypiper pipeline.

```
-h, --help            show this help message and exit
-R, --recover         Recover mode, overwrite locks
-F, --fresh-start     Fresh start mode, overwrite all
-D, --dirty           Make all cleanups manual
```

## Looper arguments
These options aren't used by Pypiper directly, but are merely provided here as a convenient way to standardize the interface for `looper.py`. If your pipeline accepts these arguments, it will be fully compatible with all the benefits of using Looper. Add these arguments with `add_pypiper_args(looper_args=True)`

```
-C CONFIG_FILE, --config CONFIG_FILE
											pipeline config file in YAML format; relative paths
											are considered relative to the pipeline script.
											defaults to rrbs.yaml
-O PARENT_OUTPUT_FOLDER, --output_parent PARENT_OUTPUT_FOLDER
											parent output directory of the project (required). The
											sample_name argument will be appended to this folder
											for output
-P NUMBER_OF_CORES, --cores NUMBER_OF_CORES
											number of cores to use for parallel processes
```

## Common arguments
Commonly used arguments that aren't directly used by Pypiper (or Looper), but which you may want to use just to keep things standard. `add_pypiper_args(common_args=True)`
```
-I INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
											one or more input files (required)
-S SAMPLE_NAME, --sample_name SAMPLE_NAME
											unique name for output subfolder and files (required)
```

## NGS arguments
Commonly used arguments that aren't directly used by Pypiper (or Looper), but which you may want to use just to keep things standard.
```
-G GENOME_ASSEMBLY, --genome GENOME_ASSEMBLY
											identifier for genome assempbly (required)
-Q SINGLE_OR_PAIRED, --single_or_paired SINGLE_OR_PAIRED
											single or paired end? default: single
```
