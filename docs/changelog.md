# Changelog

- **v0.10.0** (*Unreleased*):

    - Fixed a bug that raised exception with empty commands

- **v0.9.4** (*2019-01-31*):

    - Point release to PyPI for README rendering.

- **v0.9.3** (*2019-01-31*):

    - Simple point release update to fix PyPI landing page.

- **v0.9.2** (*2019-01-30*):

    - Never echo protected-looking attribute request.

- **v0.9.1** (*2019-01-29*):

    - Fixed a bug in NGSTk that caused errors for read counting functions on 
    MACOS. MACOS ``wc`` returns leading whitespace, which caused these functions
    to fail.

- **v0.9.0** (*2018-11-19*):

    - Use ``psutil`` to track aggregate memory usage for processes that spawn
    children. This results in accurate memory records for these processes.

    - Individual commands in a string of commands connected by shell pipes are
    now treated as individual commands, and and monitored individually for
    time and memory, and if a single component, fails, the entire string will
    fail. Previously, only the final return command was recorded, as in ``bash``.

    - Various other small improvements (like waiting checking for dynamic recover
    flags)


- **v0.8.1** (*2018-09-20*):

    - Fixed a bug that caused a problem for some pipelines adding groups of pypiper args.
    
    - Improved the `run` waiting method to immediately stop upon job
      completion, rather than minute-increment polling. This should improve
      performance particularly in pipelines with many, medium-runtime steps, and
      improve accuracy of timing profiles.


- **v0.8.0** (*2018-06-15*):

    - Implemented 'new start' mode.

    - Improved error messages and exception handling for missing child software.

    - Clarified the built-in required vs. optional args by allowing pipeline authors to specify which of the pypiper args are required. The command-line help UI now displays these correctly as 'required arguments' instead of incorrectly as 'optional arguments'.

    - Corrected the sort order of added arguments, so they are listed in the help menu more naturally.

    - Fixed a bug that caused an erroneous error message indicating missing pypiper args.

    - Clarified the license is BSD2

    - Fixed a bug that neglected to list pyyaml as a dependency

- **v0.7.2** (*2018-06-05*):

    - Implemented the 'report object' function.

    - Cleanup files are now relative, so a moved folder could still be cleaned.

    - Fixed a bug that prevented install if pypandoc was not installed

    - Fixed a bug that caused an error in containers where /proc wasn't accessible


- **v0.7.1** (*2018-02-27*):

    - Package cleanup for Pypi.

- **v0.7.0** (*2017-12-12*):

    - Standardize `NGSTk` function naming.

    - Introduce `Stage` as a model for a logically related set of pipeline processing steps.

    - Introduce `Pipeline` framework for automated processing phase execution and checkpointing.

    - Add ability to start and/or stop a pipeline at arbitrary checkpoints.

    - Introduce new state for a paused/halted pipeline.

    - Improve spawned process shutdown to avoid zombie processes.

- **v0.6** (*2017-08-24*):

    - Adds 'dynamic recovery' capability. For jobs that are terminated by an interrupt, such as a SIGINT or SIGTERM (as opposed to a failed command), pypiper will now set a dynamic recovery flags. These jobs, when restarted, will automatically pick up where they left off, without requiring any user intervention. Previously, the user would have to specify recover mode (`-R`). Now, recover mode forces a recover regardless of failure type, but interrupted pipelines will auto-recover.

    - Pypiper now appropriately adds cleanup files intermediate files for failed runs. It adds them to the cleanup script.

    - Improves error messages so only a single exception is raised with a more direct relevance to the user/

    - Pypiper will automatically remove existing flags when the run starts, eliminating the earlier issue of confusion due to multiple flags present on runs that were restarted.

    - Fixes a bug that caused a pipeline to continue if a SIGTERM is given during a process that was marked `nofail`.

    - Pypiper now can handle multiple SIGTERMs without one canceling the shutdown procedure begun by the other.

    - Major improvements to documentation and tutorials.

    - Adds `report_figure` function.

- **v0.5** (*2017-07-21*):

    - Adds preliminary support for handling docker containers

    - Updates docs, adds Hello World example

    - Adds 'waiting' flag

    - Eliminates extra spaces in reported results

    - Pypiper module is version aware

    - Updates Success time format to eliminate space

    - Improves efficiency in some ngstk merging functions

- **v0.4** (*2017-01-23*):

    - First major public release!

    - Revamps pypiper args

    - Adds parallel compression/decompression with pigz

    - Various small bug fixes and speed improvements
