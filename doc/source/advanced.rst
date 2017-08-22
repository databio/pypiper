Advanced options
=========================


Python process types: shell vs direct
****************************************************
Since Pypiper runs all your commands from within python (using the `subprocess` python module), it's nice to be aware of the two types of processes that `subprocess` allows: **direct processes** and **shell processes**.

By default, Pypiper will try to guess what kind of process you want, so for most pipelines, it's probably not necessary to understand the details in this section. However, how you write your commands has some implications for memory tracking, and advanced pipeline authors may want to control the process types that Pypiper uses, so this section covers how these subprocesses work.

**Direct process**: A direct process is one that Python executes directly, from within python. Python retains control over the process completely. Wherever possible, you should use a direct subprocess -- this has the advantage of enabling Python to monitor the memory use of the subprocess, because Python retains control over it. This the preferable way of running subprocesses in Python. The disadvantage of direct subprocesses is that you may not use shell-specific operators in a direct subprocess. For instance, if you use an asterisk (``*``) for wildcard expansion, or a bracket (``>``) for output redirection, or a pipe (``|``) to link processes -- these are commands understood by a shell like Bash, and thus, cannot be run as direct subprocesses in Python.

**Shell process**: In a shell process, Python first spawns a shell, and then runs the command in that shell. The spawned shell is then controlled by Python, but processes done by the shell are not. This allows you to use shell operators (``*``, ``|``, ``>``), but at the cost of the ability to monitor memory high water mark, because Python does not have direct control over subprocesses run inside a subshell. 

You **must** use a shell process if you are using shell operators in your command.  You can force Pypiper to use one or the other by specifying ``shell=True`` or ``shell=False`` to the ``run`` function. By default Pypiper will try to guess: if your command contains any of the shell process characters (``*``, ``|``, or ``>``), it will be run in a shell. Otherwise, it will be run as a direct subprocess. So for most purposes, you do not need to worry about this at all, but you may want to write your commands to minimize shell operators if you are interested in the memory monitoring features of Pypiper.

