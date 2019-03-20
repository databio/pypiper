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


### atexit\_register
Convenience alias to register exit functions without having to import atexit in the pipeline.
```python
def atexit_register(self, *args):
```




### callprint
Prints the command, and then executes it, then prints the memory use and return code of the command.

Uses python's subprocess.Popen() to execute the given command. The shell argument is simply
passed along to Popen(). You should use shell=False (default) where possible, because this enables memory
profiling. You should use shell=True if you require shell functions like redirects (>) or stars (*), but this
will prevent the script from monitoring memory use. The pipes (|) will be used to split the command into
subprocesses run within python, so the memory profiling is possible.
cmd can also be a series (a dict object) of multiple commands, which will be run in succession.
```python
def callprint(self, cmd, shell=None, lock_file=None, nofail=False, container=None):
```

**Parameters:**

- `cmd` -- `str | Iterable[str]`:  Bash command(s) to be run.
- `lock_file` -- `str`:  a lock file name
- `nofail` -- `bool`:  FalseNofail can be used to implement non-essential parts of the pipeline; if these processes fail, they will not cause the pipeline to bail out.
- `shell` -- `bool`:  None (will tryto determine based on the command)
- `container` -- ``:  Named Docker container in which to execute.
- `container` -- ``:  str




### checkprint
Just like callprint, but checks output -- so you can get a variable in python corresponding to the return value of the command you call. This is equivalent to running subprocess.check_output() instead of subprocess.call().
```python
def checkprint(self, cmd, shell=None, nofail=False, errmsg=None):
```

**Parameters:**

- `cmd` -- `str | Iterable[str]`:  Bash command(s) to be run.
- `shell` -- `bool | str`:  If command requires should be run in its own shell. Optional.Default: "guess" -- `run()` will try to guess if the command should be run in a shell (based on the presence of a pipe (|) or redirect (>), To force a process to run as a direct subprocess, set `shell` to False; to force a shell, set True.
- `nofail` -- `bool`:  FalseNofail can be used to implement non-essential parts of the pipeline; if these processes fail, they will not cause the pipeline to bail out.
- `errmsg` -- `str`:  Message to print if there's an error.




### clean\_add
Add files (or regexs) to a cleanup list, to delete when this pipeline completes successfully. When making a call with run that produces intermediate files that should be deleted after the pipeline completes, you flag these files for deletion with this command. Files added with clean_add will only be deleted upon success of the pipeline.
```python
def clean_add(self, regex, conditional=False, manual=False):
```

**Parameters:**

- `regex` -- `str`:   A unix-style regular expression that matches files to delete(can also be a file name).
- `conditional` -- `bool`:  True means the files will only be deleted if no otherpipelines are currently running; otherwise they are added to a manual cleanup script called {pipeline_name}_cleanup.sh
- `manual` -- `bool`:  True means the files will just be added to a manual cleanup script.




### complete
Stop a completely finished pipeline.
```python
def complete(self):
```




### completed
Is the managed pipeline in a completed state?
```python
def completed:
```

**Returns:**

`bool`:  Whether the managed pipeline is in a completed state.




### fail\_pipeline
If the pipeline does not complete, this function will stop the pipeline gracefully. It sets the status flag to failed and skips the normal success completion procedure.
```python
def fail_pipeline(self, e, dynamic_recover=False):
```

**Parameters:**

- `e` -- `Exception`:  Exception to raise.
- `dynamic_recover` -- `bool`:  Whether to recover e.g. for job termination.




### failed
Is the managed pipeline in a failed state?
```python
def failed:
```

**Returns:**

`bool`:  Whether the managed pipeline is in a failed state.




### flag\_file\_path
Create path to flag file based on indicated or current status.

Internal variables used are the pipeline name and the designated
pipeline output folder path.
```python
def flag_file_path(self, status=None):
```

**Parameters:**

- `status` -- `str`:  flag file type to create, default to current status


**Returns:**

`str`:  path to flag file of indicated or current status.




### get\_container
```python
def get_container(self, image, mounts):
```



### get\_stat
Returns a stat that was previously reported. This is necessary for reporting new stats that are derived from two stats, one of which may have been reported by an earlier run. For example, if you first use report_result to report (number of trimmed reads), and then in a later stage want to report alignment rate, then this second stat (alignment rate) will require knowing the first stat (number of trimmed reads); however, that may not have been calculated in the current pipeline run, so we must retrieve it from the stats.tsv output file. This command will retrieve such previously reported stats if they were not already calculated in the current pipeline run.
```python
def get_stat(self, key):
```

**Parameters:**

- `key` -- ``:  key of stat to retrieve




### halt
Stop the pipeline before completion point.
```python
def halt(self, checkpoint=None, finished=False, raise_error=True):
```

**Parameters:**

- `checkpoint` -- `str`:  Name of stage just reached or just completed.
- `finished` -- `bool`:  Whether the indicated stage was just finished(True), or just reached (False)
- `raise_error` -- `bool`:  Whether to raise an exception to trulyhalt execution.




### halted
Is the managed pipeline in a paused/halted state?
```python
def halted:
```

**Returns:**

`bool`:  Whether the managed pipeline is in a paused/halted state.




### has\_exit\_status
Has the managed pipeline been safely stopped?
```python
def has_exit_status:
```

**Returns:**

`bool`:  Whether the managed pipeline's status indicates that ithas been safely stopped.




### is\_running
Is the managed pipeline running?
```python
def is_running:
```

**Returns:**

`bool`:  Whether the managed pipeline is running.




### make\_sure\_path\_exists
Creates all directories in a path if it does not exist.
```python
def make_sure_path_exists(self, path):
```

**Parameters:**

- `path` -- `str`:  Path to create.


**Raises:**

- `Exception`:  if the path creation attempt hits an error witha code indicating a cause other than pre-existence.




### remove\_container
```python
def remove_container(self, container):
```



### report\_figure
Writes a string to self.pipeline_figures_file.
```python
def report_figure(self, key, filename, annotation=None):
```

**Parameters:**

- `key` -- `str`:  name (key) of the figure
- `filename` -- `str`:  relative path to the file (relative to parent output dir)
- `annotation` -- `str`:  By default, the figures will be annotated with the pipelinename, so you can tell which pipeline records which figures. If you want, you can change this.




### report\_object
Writes a string to self.pipeline_objects_file. Replacement of report_figure
```python
def report_object(self, key, filename, anchor_text=None, anchor_image=None, annotation=None):
```

**Parameters:**

- `key` -- `str`:  name (key) of the object
- `filename` -- `str`:  relative path to the file (relative to parent output dir)
- `anchor_text` -- `str`:  text used as the link anchor test or caption torefer to the object. If not provided, defaults to the key.
- `anchor_image` -- `str`:  a path to an HTML-displayable image thumbnail (so,.png or .jpg, for example). If a path, the path should be relative to the parent output dir.
- `annotation` -- `str`:  By default, the figures will be annotated with thepipeline name, so you can tell which pipeline records which figures. If you want, you can change this.




### report\_result
Writes a string to self.pipeline_stats_file.
```python
def report_result(self, key, value, annotation=None):
```

**Parameters:**

- `key` -- `str`:  name (key) of the stat
- `annotation` -- `str`:  By default, the stats will be annotated with the pipelinename, so you can tell which pipeline records which stats. If you want, you can change this; use annotation='shared' if you need the stat to be used by another pipeline (using get_stat()).




### run
The primary workhorse function of PipelineManager, this runs a command.

This is the command  execution function, which enforces
race-free file-locking, enables restartability, and multiple pipelines
can produce/use the same files. The function will wait for the file
lock if it exists, and not produce new output (by default) if the
target output file already exists. If the output is to be created,
it will first create a lock file to prevent other calls to run
(for example, in parallel pipelines) from touching the file while it
is being created. It also records the memory of the process and
provides some logging output.
```python
def run(self, cmd, target=None, lock_name=None, shell=None, nofail=False, clean=False, follow=None, container=None):
```

**Parameters:**

- `cmd` -- `str | list[str]`:  Shell command(s) to be run.
- `target` -- `str | Sequence[str]`:  Output file(s) to produce, optional
- `lock_name` -- `str`:  Name of lock file. Optional.
- `shell` -- `bool`:  If command requires should be run in its own shell.Optional. Default: None --will try to determine whether the command requires a shell.
- `nofail` -- `bool`:  Whether the pipeline proceed past a nonzero return froma process, default False; nofail can be used to implement non-essential parts of the pipeline; if a 'nofail' command fails, the pipeline is free to continue execution.
- `errmsg` -- `str`:  Message to print if there's an error.
- `clean` -- `bool`:  True means the target file will be automatically addedto an auto cleanup list. Optional.
- `follow` -- `callable`:  Function to call after executing (each) command.
- `container` -- `str`:  Name for Docker container in which to run commands.


**Returns:**

`int`:  Return code of process. If a list of commands is passed,this is the maximum of all return codes for all commands.




### set\_status\_flag
Configure state and files on disk to match current processing status.
```python
def set_status_flag(self, status):
```

**Parameters:**

- `status` -- `str`:  Name of new status designation for pipeline.




### start\_pipeline
Initialize setup. Do some setup, like tee output, print some diagnostics, create temp files. You provide only the output directory (used for pipeline stats, log, and status flag files).
```python
def start_pipeline(self, args=None, multi=False):
```




### stop\_pipeline
Terminate the pipeline.

This is the "healthy" pipeline completion function.
The normal pipeline completion function, to be run by the pipeline
at the end of the script. It sets status flag to completed and records 
some time and memory statistics to the log file.
```python
def stop_pipeline(self, status='completed'):
```




### time\_elapsed
Returns the number of seconds that have elapsed since the time_since parameter.
```python
def time_elapsed(self, time_since):
```

**Parameters:**

- `time_since` -- `float`:  Time as a float given by time.time().




### timestamp
Print message, time, and time elapsed, perhaps creating checkpoint.

This prints your given message, along with the current time, and time
elapsed since the previous timestamp() call.  If you specify a
HEADING by beginning the message with "###", it surrounds the message
with newlines for easier readability in the log file. If a checkpoint
is designated, an empty file is created corresponding to the name
given. Depending on how this manager's been configured, the value of
the checkpoint, and whether this timestamp indicates initiation or
completion of a group of pipeline steps, this call may stop the
pipeline's execution.
```python
def timestamp(self, message='', checkpoint=None, finished=False, raise_error=True):
```

**Parameters:**

- `message` -- `str`:  Message to timestamp.
- `checkpoint` -- `str`:  Name of checkpoint; this tends to be somethingthat reflects the processing logic about to be or having just been completed. Provision of an argument to this parameter means that a checkpoint file will be created, facilitating arbitrary starting and stopping point for the pipeline as desired.
- `finished` -- `bool`:  Whether this call represents the completion of aconceptual unit of a pipeline's processing
- `raise_error` -- ``:  Whether to raise exception ifcheckpoint or current state indicates that a halt should occur.



