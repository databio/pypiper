

The ``follow`` argument
================================================================================

.. _the_follow_argument:

The ``PipelineManager.run`` function has an optional argument named ``follow`` that is useful for checking or reporting results from a command. To the ``follow`` argument you must pass a python function (which may be either a defined function or a ``lambda`` function). These *follow functions* are then coupled to the command that is run; the follow function will be called by python **if and only if** the command is run. 

Why is this useful? The major use cases are QC checks and reporting results. We use a folllow function to  run a QC check to make sure processes did what we expect, and then to report that result to the ``stats`` file. We only need to check the result and report the statistic once, so it's best to put these kind of checks in a ``follow`` function. Often, you'd like to run a function to examine the result of a command, but you only want to run that once, *right after the command that produced the result*. For example, counting the number of lines in a file after producing it, or counting the number of reads that aligned right after an alignment step. You want the counting process coupled to the alignment process, and don't need to re-run the counting every time you restart the pipeline. Because pypiper is smart, it will not re-run the alignment once it has been run; so there is no need to re-count the result on every pipeline run! 

*Follow functions* let you avoid running unnecessary processes repeatedly in the event that you restart your pipeline multiple times (for instance, while debugging later steps in the pipeline).
