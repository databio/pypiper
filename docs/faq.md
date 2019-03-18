# FAQ

- **How can I run my pipeline on more than 1 sample?**
	Pypiper only handles individual-sample pipelines. To run it on multiple samples, write a loop, or use `Looper <http://looper.readthedocs.io/>`_. Dividing multi-sample handling from individual sample handling is a conceptual advantage that allows us to write a nice, universal, generic sample-handler that you only have to learn once.

- **What cluster resources can pypiper use?** 
	Pypiper is compute-agnostic. You run it wherever you want; If you want a nice way to submit pipelines for samples any cluster manager, check out `Looper <http://looper.readthedocs.io/>`_.

- **What does it mean for a sample to be in the "waiting" state?**
	Waiting means it encountered a file lock, but no recovery flag. So the pipeline thinks a process (from another run or another process) is currently writing that file. It periodically checks for the lock file to disappear, and assumes that the other process will unlock the file when finished. If the lock was left by a previous failed run, then it will just wait forever. This is what recover mode (``-R``) is intended for, if pipelines fail.

- **What is the 'elapsed time' in output?**
	The "elapsed" time is referring to the amount of time since the preceding timestamp, not since the start of the pipeline. Timestamps are all displayed with a flag: ``_TIME_``. The total cumulative time for the pipeline is displayed only at the end.

- **How should I run a QC step to check results of one of my commands?**
	Usually, you only want to run a QC step if the result was created in the same pipeline run. There's no need to re-run that step if you have to restart the pipeline due to an error later on. If you use ``run()`` for these steps, then they'll need to run each time the pipeline runs. Instead, this is exactly why we created :ref:`the follow argument <the_follow_argument>`. This option lets you couple a QC step to a ``run()`` call, so it only gets excecuted when it is required.