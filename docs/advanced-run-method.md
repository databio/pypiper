# Avanced run method options

The `PipelineManager.run()` function is the core part of `pypiper`. In its simplest case, all you need to provide is a command to run. But it can be much more powerful than that if you make use of these additional arguments:

## The `cmd` argument

Normally you just pass a string, but you can also pass a list of commands to `run`, like this:

```
pm.run([cmd1, cmd2, cmd3])
```

Pypiper will treat these commands as a group, running each one in turn (and monitoring them individually for time and memory use). The difference in doing it this way, rather than 3 separate calls to `run()` is that if the series does not complete, the entire series will be re-run. This is therefore useful to piece together commands that must all be run together.

## The `target` and `lock_name` arguments

If you provide a `target` file, then `pypiper` will first check to see if that target exists, and only run the `command` if the `target` does not exist. To prevent two pipelines from running commands on the same target, `pypiper` will automatically derive a lock file name from your target file. You can use the `lock_name` argument to override this default. If you do not provide a `target`, then you will need to provide a `lock_name` argument because `pypiper` will not be able to derive one automatically.

## The `shell` argument: Python subprocess types

Since Pypiper runs all your commands from within python (using the `subprocess` python module), it's nice to be aware of the two types of processes that `subprocess` allows: **direct processes** and **shell processes**.

By default, Pypiper will try to guess what kind of process you want, so for most pipelines, it's probably not necessary to understand the details in this section. However, how you write your commands has some implications for memory tracking, and advanced pipeline authors may want to control the process types that Pypiper uses, so this section covers how these subprocesses work.

**Direct process**: A direct process is one that Python executes directly, from within python. Python retains control over the process completely. Wherever possible, you should use a direct subprocess -- this has the advantage of enabling Python to monitor the memory use of the subprocess, because Python retains control over it. This the preferable way of running subprocesses in Python. The disadvantage of direct subprocesses is that you may not use shell-specific operators in a direct subprocess. For instance, if you use an asterisk (`*`) for wildcard expansion, or a bracket (`>`) for output redirection, or a pipe (`|`) to link processes -- these are commands understood by a shell like Bash, and thus, cannot be run as direct subprocesses in Python.

**Shell process**: In a shell process, Python first spawns a shell, and then runs the command in that shell. The spawned shell is then controlled by Python, but processes done by the shell are not. This allows you to use shell operators (`*`, `|`, `>`), but at the cost of the ability to monitor memory high water mark, because Python does not have direct control over subprocesses run inside a subshell. 

You **must** use a shell process if you are using shell operators in your command.  You can force Pypiper to use one or the other by specifying `shell=True` or `shell=False` to the `run` function. By default Pypiper will try to guess: if your command contains any of the shell process characters (`*`, `|`, or `>`), it will be run in a shell. Otherwise, it will be run as a direct subprocess. So for most purposes, you do not need to worry about this at all, but you may want to write your commands to minimize shell operators if you are interested in the memory monitoring features of Pypiper.

## The `nofail` argument

By default, a command that fails will cause the entire pipeline to halt. If you want to provide an optional command that *should not* halt the pipeline upon failure, set `nofail=True`. `nofail` can be used to implement non-essential parts of the pipeline.
 


## The `follow` argument

The `PipelineManager.run` function has an optional argument named `follow` that is useful for checking or reporting results from a command. To the `follow` argument you must pass a python function (which may be either a defined function or a `lambda` function). These *follow functions* are then coupled to the command that is run; the follow function will be called by python **if and only if** the command is run. 

Why is this useful? The major use cases are QC checks and reporting results. We use a folllow function to  run a QC check to make sure processes did what we expect, and then to report that result to the `stats` file. We only need to check the result and report the statistic once, so it's best to put these kind of checks in a `follow` function. Often, you'd like to run a function to examine the result of a command, but you only want to run that once, *right after the command that produced the result*. For example, counting the number of lines in a file after producing it, or counting the number of reads that aligned right after an alignment step. You want the counting process coupled to the alignment process, and don't need to re-run the counting every time you restart the pipeline. Because pypiper is smart, it will not re-run the alignment once it has been run; so there is no need to re-count the result on every pipeline run! 

*Follow functions* let you avoid running unnecessary processes repeatedly in the event that you restart your pipeline multiple times (for instance, while debugging later steps in the pipeline).

## The `container` argument

If you specify a string here, `pypiper` will wrap the command in a `docker run` call using the given `container` image name.
