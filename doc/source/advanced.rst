Advanced
=========================


Python process types: Shell vs direct
*************
Since Pypiper runs all your commands from within python (using the `subprocess` python module), you need to be aware of the two types of processes that `subprocess` can handle: Shell process and direct processes. Choosing between the two is as simple as specifying `shell=True` (the default is a direct subprocess) to `run()`, which is directly passed on as the `shell` parameter to `subprocess.Popen()`. The difference is that Python runs the process either as a subprocess directly from within python, or it first spawns a shell, and then runs the subprocess in the shell. When should you use each?

**Direct process**: For most use cases, you should simply use a direct subprocess (the default) -- this has the advantage of enabling Python to monitor the memory use of the subprocess, because Python retains control over it. This is considered by the community to be the preferable way of running subprocesses in Python.

**Shell process**: You must use a `shell=True` process if you are using shell operators in your command. For instance, if you use an asterisk (`*`) for wildcard expansion, or a bracket (`>`) for output redirection, or a pipe (`|`) to link processes -- these are commands understood by a shell like Bash, and thus, cannot be run as direct subprocesses by Python, but must be run in a shell. `Subprocess` (and therefore `Pypiper`) can understand these commands just fine, but you lose the ability to monitor memory high water mark because Python does not have direct control over subprocesses run inside a subshell.

