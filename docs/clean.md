# Cleaning up intermediate files

Many pipelines produce intermediate files along the way. Should you retain these files or delete them?

On the one hand, you may not necessarily want to delete them *immediately* after creating them, because what if a later pipeline step fails and you need to inspect an intermediate file? On the other hand, you may not want those intermediate files sticking around forever because they waste valuable disk space.

Pypiper solves this problem with the concept of a *clean list*. The clean list is simply a list of files that are flagged for eventual cleanup. A pipeline developer adds to this list using `pm.clean_add(filename)`. Files on the clean list are *not* cleaned immediately; instead, they are **removed as soon as the pipeline is completed successfully** (in other words, after `pm.complete_pipeline()` is called). The advantage is that intermediate files will always be available as long as a pipeline has not completed successfully.

In case a user of a pipeline instead wants to retain these files indefinitely, he or she may simply add `--dirty` when invoking the pipeline script. This instructs pypiper to *not* clean the intermediate files, even after a successful pipeline run. In this case, `pypiper` will produce a shell script (`clean.sh`), which can be run to remove all flagged files at a later point.
