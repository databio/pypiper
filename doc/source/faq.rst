
FAQ
=========================

- **How can I run my pipeline on more than 1 sample?**
	Pypiper only handles individual-sample pipelines. To run it on multiple samples, write a loop, or use `Looper <http://looper.readthedocs.io/>`_. Dividing multi-sample handling from individual sample handling is a conceptual advantage that allows us to write a nice, universal, generic sample-handler that you only have to learn once.

- **What cluster resources can pypiper use?** 
	PyPiper is compute-agnostic. You run it wherever you want; If you want a nice way to submit pipelines for samples any cluster manager, check out `Looper <http://looper.readthedocs.io/>`_.

- **What does it mean for a sample to be in the "waiting" state?**
	Waiting means it encountered a file lock, but no recovery flag. So the pipeline thinks a process (from another run or another process) is currently writing that file. It periodically checks for the lock file to disappear, and assumes that the other process will unlock the file when finished. If the lock was left by a previous failed run, then it will just wait forever. This is what recover mode (``-R``) is intended for, if pipelines fail.

- **What is the 'elapsed time' in output?**
	The "elapsed" time is referring to the amount of time since the preceding timestamp, not since the start of the pipeline. Timestamps are all displayed with a flag: ``_TIME_``. The total cumulative time for the pipeline is displayed only at the end.
