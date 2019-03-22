# Package pypiper Documentation

## Class PipelineManager
Base class for instantiating a PipelineManager object, the main class of Pypiper.

**Parameters:**

- `name` -- `str`:  Choose a name for your pipeline;it's used to name the output files, flags, etc.
- `outfolder` -- `str`:  Folder in which to store the results.
- `args` -- `argparse.Namespace`:  Optional args object from ArgumentParser;Pypiper will simply record these arguments from your script
- `multi` -- `bool`:  Enables running multiple pipelines in one scriptor for interactive use. It simply disables the tee of the output, so you won't get output logged to a file.
- `dirty` -- `bool`:  Overrides the pipeline's clean_add()manual parameters, to *never* clean up intermediate files automatically. Useful for debugging; all cleanup files are added to manual cleanup script.
- `recover` -- `bool`:  Specify recover mode, to overwrite lock files.If pypiper encounters a locked target, it will ignore the lock and recompute this step. Useful to restart a failed pipeline.
- `new_start` -- `bool`:  start over and run every command even if output exists
- `force_follow` -- `bool`:  Force run all follow functionseven if  the preceding command is not run. By default, following functions  are only run if the preceding command is run.
- `cores` -- `int`:  number of processors to use, default 1
- `mem` -- `str`:  amount of memory to use. Default units are megabytes unlessspecified using the suffix [K|M|G|T]."
- `config_file` -- `str`:  path to pipeline configuration file, optional
- `output_parent` -- `str`:  path to folder in which output folder will live
- `overwrite_checkpoints` -- `bool`:  Whether to override the stage-skippinglogic provided by the checkpointing system. This is useful if the calls to this manager's run() method will be coming from a class that implements pypiper.Pipeline, as such a class will handle checkpointing logic automatically, and will set this to True to protect from a case in which a restart begins upstream of a stage for which a checkpoint file already exists, but that depends on the upstream stage and thus should be rerun if it's "parent" is rerun.


**Raises:**

- `TypeError`:  if start or stop point(s) are provided both directly andvia args namespace, or if both stopping types (exclusive/prospective and inclusive/retrospective) are provided.


### halted
Is the managed pipeline in a paused/halted state?
```python
def halted(self)
```

**Returns:**

`bool`:  Whether the managed pipeline is in a paused/halted state.



