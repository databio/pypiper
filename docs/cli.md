# Command-line arguments

Your final pypiper pipeline will be a python script that a pipeline user will invoke on the command-line. You will likely need to allow the user to change some parameters on the command line, and to take full advantage of Pypiper (make your pipeline recoverable, etc.), you will need to add command-line options to your pipeline that change pypiper's settings as well. Pypiper uses the typical Python [argparse module](https://docs.python.org/2/library/argparse.html) to define command-line arguments to your pipeline, and offers a series of built-in functions to help you populate your pipeline's `ArgumentParser` with pypiper-specific options.

You can use an ArgumentParser as usual, adding whatever arguments you like. Then, you add Pypiper args to your parser with the function `add_pypiper_args()`, and pass command-line options and arguments to your `PipelineManager`, like this:

```{python}
import pypiper, os, argparse
parser = ArgumentParser(description='Write a short description here')

# add any custom args here
# e.g. parser.add_argument('--foo', help='foo help')

# once you've established all your custom arguments, we can add the default
# pypiper arguments to your parser like this:

parser = pypiper.add_pypiper_args(parser)

# Then, pass the args parsed along to the PipelineManger

args = parser.parse_args()

pipeline = pypiper.PipelineManager(name="my_pipeline", outfolder="out", \
                    args=args)
```

Once you've added pypiper arguments, your pipeline will then enable a few built-in arguments: `--recover`, `--follow`, and `--dirty`, for example. As a side bonus, all arguments (including any of your custom arguments) will be recorded in the log outputs. 

That's the basics. But you can customize things for more efficiency using a simple set of pre-built args and groups of args in pypiper:


# Universal pypiper options

With that said, there are a few universal (Pypiper-added) options that are frequently (but not necessarily always) honored by pypiper pipelines. These default pypiper arguments are detailed below:

  - `-R, --recover`
    Recover mode, overwrite locks. This argument will tell pypiper to recover from a failed previous run. Pypiper will execute commands until it encounters a locked file, at which point it will re-execute the failed command and continue from there.

  - `-F, --follow`
    Force run follow-functions. By default, follow-functions are only run if their corresponding `run` command was run; with this option you can force all follow functions to run. This is useful for regenerating QC data on existing output. For more details, see :ref:`the follow argument <the_follow_argument>`.

  - `-D, --dirty`
    Make all cleanups manual. By default, pypiper pipelines will delete any intermediate files. For debugging, you may want to turn this option off -- you can do that by specifying **dirty mode**.

  - `-N, --new-start`
    New start mode. This flag will tell pypiper to start over, and run every command, even if its target output already exists.


## Customizing `add_pypiper_args()`


There are two ways to modulate the arguments added by `add_pypiper_args()` function: the `groups` argument, which lets you add argument groups; or the `args` argument, which lets you add arguments individually. By default, `add_pypiper_args()` add all arguments listed in the `pypiper` group. You may instead pass a list of one or more of these groups of arguments (to `groups`) or individual arguments (to `args`) to customize exactly the set of built-in options your pipeline implements.

For example, `parser.add_pypiper_args(parser, groups=['pypiper', 'common'])` will add all arguments listed under `pypiper` and `common` below:


## Built-in arguments accessed with `add_pypiper_args()`

Individual arguments that are understood and used by pypiper:

- `-R, --recover`: for a failed pipeline run, start off at the last successful step. 
- `-N, --new-start`: Just recreate everything, even if it exists.
- `-D, --dirty`: Disables automatic cleaning of temporary files, so all intermediate files will still exist after a pipeline run (either successful or failed). Useful for debugging a pipeline even if it succeeds.
- `-F, --follow`: Runs all `follow-functions`, regardless of whether the accompanying command is run.
- `-C, --config`: Pypiper pipeline config yaml file.

Individual arguments just provided for convenience and standardization:
- `-S, --sample-name`: name of the sample
- `-I, --input`: primary input file (e.g. read1)
- `-I2, --input2`: secondary input file (e.g. read2)
- `-O, --output-parent`: parent folder for pipeline results (the pipeline will use this as the parent directory for a folder named `sample-name`)
- `-P, --cores`: Number of cores to use
- `-M, --mem`: Amount of memory in megabytes
- `-G, --genome`: Reference genome assembly (e.g. `hg38`)
- `-Q, --simple-or-paired`: For sequencing data, is input single-end or paired-end?

## Pre-built collections of arguments added via `groups`:

- pypiper: `recover`, `new-start`, `dirty`, `follow`
- common: `input`, `sample-name`
- config: `config`
- resource: `mem`, `cores`
- looper: `config`, `output-parent`, `mem`, `cores`
- ngs: `input`, `sample-name`, `input2`, `genome`, `single-or-paired`


## Specifying required built-in arguments

If you're using the built-in arguments, you may want to module which are required and which are not. That way, you can piggyback on how `ArgumentParser` handles required arguments very nicely -- if the user does not specify a required argument, the pipeline will automatically prompt with usage instructions.

By default, built-in arguments are not flagged as required, but you can pass a list of required built-ins to the `required` parameter, like `add_pypiper_args(parser, args=["sample-name"], required=["sample-name"])`.


## Examples

    import pypiper, os, argparse
    parser = ArgumentParser(description='Write a short description here')

    # add just arguments from group `pypiper`
    parser = pypiper.add_pypiper_args(parser, groups=["pypiper"])

    # add just arguments from group `common`
    parser = pypiper.add_pypiper_args(parser, groups=["common"])    

    # add arguments from two groups
    parser = pypiper.add_pypiper_args(parser, groups=["common", "resources"],
                                        required=["sample-name", "output-parent"])

    # add individual argument
    parser = pypiper.add_pypiper_args(parser, args=["genome"])

    # add some groups and some individual arguments
    parser = pypiper.add_pypiper_args(parser, args=["genome"], groups=["looper", "ngs"])
