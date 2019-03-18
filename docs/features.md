# Pypiper features at-a-glance

![](img/simplicity.svg) **Simplicity**

Pipelines are simple both to use and to develop. A pypiper pipeline is nothing more than a python script. You run it on the command line like you would any other python script. The basic documentation is just a few pages. It should only take you 15 minutes to write your first pipeline. 

![](img/restartability.svg) **Restartability**

Commands check for their targets and only run if the target needs to be created. This provides computational advantages, and also means the pipeline will pick up where it left off in case it needs to be restarted or extended.

![](img/protection.svg) **File integrity protection**

Pypiper uses automatic file locks. This ensures that tasks complete, and pipelines never continue with half-finished analysis. It also ensures that multiple pipeline runs will not interfere with one another -even if the steps are identical and produce the same files.

![](img/logging.svg) **Copious logging**

Pypiper automatically prints output to screen and also stores it in a log file, so all subprocess output is captured permanently. It also provides copious information on versions, compute host, and easy timestamping.

![](img/memory.svg) **Memory use monitoring**

Processes are polled for memory use, allowing you to more accurately gauge your future memory requirements.

![](img/job_status.svg) **Job status monitoring**

Pypiper automatically creates status flag files, so you can summarize the current state (`running`, `failed`, or `completed`) of hundreds of jobs simultaneously.

![](img/reports.svg) **Easy result reports**

Pypiper provides functions to put key-value pairs into an easy-to-parse stats file, making it easy to summarize your pipeline results.

![](img/error.svg) **Robust error handling**

Pypiper closes pipelines gracefully on interrupt or termination signals, converting the status to `failed`. By default, a process that returns a nonzero value halts the pipeline, unlike in bash, where by default the pipeline would continue using an incomplete or failed result. This behavior can be overridden as desired with a single parameter.

![](img/recovery.svg) **Dynamic recovery**

If a job is interrupted (with SIGINT or SIGTERM), either from a user or by a cluster resource manager, pypiper will set a `dynamic recovery` flag. The next time the run is started, it will automatically pick up where it left off. This makes pypiper pipelines `automatically pre-emption ready`, so they can be immediately deployed on servers where jobs may be pre-empted.
