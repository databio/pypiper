API
===

This page contains a comprehensive list of all functions within ``pypiper``.
Docstrings should provide sufficient understanding for any individual function.

pypiper.PipelineManager
-----------------------

.. currentmodule:: pypiper.pypiper.PipelineManager
.. autosummary::
   callprint
   checkprint
   clean_add
   cleanup
   create_file
   create_file_racefree
   exit_handler
   fail_pipeline
   ignore_interrupts
   kill_child_process
   make_sure_path_exists
   memory_usage
   report_command
   report_profile
   report_result
   run
   set_status_flag
   signal_int_handler
   signal_term_handler
   start_pipeline
   stop_pipeline
   time_elapsed
   timestamp
   wait_for_lock
   wait_for_process

pypiper.ngstk
-----------------------

.. currentmodule:: pypiper.ngstk.NGSTk
.. autosummary::
   addTrackToHub
   AnnotatePeaks
   bam2fastq
   bam_conversions
   bam_to_fastq
   bamToBed
   bamToBigWig
   bowtie2Map
   calculateFRiP
   centerPeaksOnMotifs
   count_fail_reads
   count_flag_reads
   count_lines
   count_mapped_reads
   count_multimapping_reads
   count_reads
   count_unique_mapped_reads
   count_unique_reads
   count_uniquelymapping_reads
   fastqc
   fastQC
   filterPeaksMappability
   filterReads
   genomeWideCoverage
   getFragmentSizes
   getFRiP
   getPeakNumber
   getReadType
   homerFindMotifs
   htSeqCount
   indexBam
   kallisto
   linkToTrackHub
   macs2CallPeaks
   macs2CallPeaksATACSeq
   macs2PlotModel
   make_sure_path_exists
   makeDir
   markDuplicates
   merge_bams
   mergeBams
   moveFile
   parseBowtieStats
   parseDuplicateStats
   parseQC
   peakTools
   picardMarkDuplicates
   plotInsertSizesFit
   preseq_coverage
   preseq_curve
   preseq_extrapolate
   removeDuplicates
   removeFile
   sam_conversions
   samtools_index
   samtools_view
   shiftReads
   skewer
   slurmFooter
   slurmHeader
   slurmSubmitJob
   sortIndexBam
   sppCallPeaks
   topHatMap
   trimmomatic
   zinbaCallPeaks


Definitions
-----------------------

.. automodule:: pypiper
   :members:

.. automodule:: pypiper.pypiper
   :members:

.. automodule:: pypiper.ngstk
   :members:
