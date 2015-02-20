#!/usr/bin/python2.7
"""
RNA BitSeq pipeline
documentation.
"""

from argparse import ArgumentParser
import pipetk
import os
import os.path
import sys
from subprocess import call
import subprocess

import ngstk

from datetime import datetime

#Set pipeline ID
pipeline = "rnaBitSeq"


# Argument Parsing
# #######################################################################################
parser = ArgumentParser(description='RNA BitSeq pipeline.')

parser.add_argument('-i', '--unmapped-bam',
                    default="/fhgfs/groups/lab_bock/jklughammer/projects/otherProjects/CORE-seq/titration/CORE/unmapped_bam/BSF_0131_C5FD6ACXX_8__CORE_K562_500_1_sub.bam",
                    nargs="+", dest='unmapped_bam', help="Input unmapped bam file(s))")

parser.add_argument('-s', '--sample-name', default="default",
                    dest='sample_name', type=str, help='Sample Name')

parser.add_argument('-r', '--project-root', default="/fhgfs/groups/lab_bock/shared/COREseq/results_pipeline/",
                    dest='project_root', type=str,
                    help='Directory in which the project will reside. Default=/fhgfs/groups/lab_bock/shared/COREseq/results_pipeline/')

parser.add_argument('-g', '--genome', default="hg19_cdna",
                    dest='genome_assembly', type=str, help='Genome Assembly')
parser.add_argument('-e', '--ercc', default="ERCC92",
                    dest='ERCC_assembly', type=str, help='ERCC Assembly')

parser.add_argument('--no-checks', dest='sanity_check', action='store_false', default=True)
parser.add_argument('-p', '--paired-end', dest='paired_end', action='store_true', default=True, help='Paired End Mode')
parser.add_argument('-q', '--single-end', dest='paired_end', action='store_false', help='Single End Mode')

args = parser.parse_args()


# Merging
########################################################################################
# If 2 unmapped bam files are given, then these are to be merged.
# Must be done here to initialize the sample name correctly
merge = False
if (len(args.unmapped_bam) > 1):
    merge = True
    if (args.sample_name == "default"):
        args.sample_name = "merged";
else:
    if (args.sample_name == "default"):
        args.sample_name = os.path.splitext(os.path.basename(args.unmapped_bam[0]))[0]

# Set up environment path variables
########################################################################################
# Set up an container class to hold paths
class Container:
    pass


paths = Container()
paths.picard_dir = "/fhgfs/groups/lab_bsf/meth-seq/src/tools/picard-tools-1.100"
paths.trimmomatic_jar = "/cm/shared/apps/trimmomatic/0.32/trimmomatic-0.32-epignome.jar"
paths.scripts_dir = os.path.dirname(os.path.realpath(__file__))
paths.adapter_file = "/fhgfs/groups/lab_bock/jklughammer/projects/COREseq/data/adapters/epignome_adapters_2_add.fa"
paths.base_folder = "/fhgfs/groups/lab_bock/jklughammer/projects/comparative_epigenomics"
paths.bowtie_indexed_genome = paths.base_folder + "/genomes/" + args.genome_assembly + "/forBowtie1/" + args.genome_assembly
paths.bowtie_indexed_ERCC = paths.base_folder + "/genomes/" + args.ERCC_assembly + "/forBowtie1/" + args.ERCC_assembly
#paths.bismark_spikein_genome = "/fhgfs/groups/lab_bsf/resources/ref_genome_5hmC"
paths.bismark_spikein_genome = "/fhgfs/groups/lab_bsf/resources/ref_genome_k1_k3"
paths.src_folder = "/fhgfs/groups/lab_bock/jklughammer/pipelines"
paths.ref_genome_fasta = paths.base_folder + "/genomes/" + args.genome_assembly + "/" + args.genome_assembly + ".fa"
paths.ref_ERCC_fasta = paths.base_folder + "/genomes/" + args.ERCC_assembly + "/" + args.ERCC_assembly + ".fa"
paths.ref_genome = paths.base_folder + "/genomes/"
paths.pipeline_outfolder = os.path.join(args.project_root + args.sample_name + "/")
paths.log_file = paths.pipeline_outfolder + "rnaBitSeq.log.md"
paths.bowtie1 = "/cm/shared/apps/bowtie/1.1.1/bin/bowtie"
paths.pipe_stats = paths.pipeline_outfolder + "/" + "stats_rnaBitSeq.txt"

# Run some initial setting and logging code to start the pipeline.
start_time = pipetk.start_pipeline(paths, args , pipeline )

print("N input bams:\t\t" + str(len(args.unmapped_bam)))
print("Sample name:\t\t" + args.sample_name)

sample_merged_bam = args.sample_name + ".merged.bam"
pipetk.make_sure_path_exists(paths.pipeline_outfolder + "unmapped_bam/")

if merge and not os.path.isfile(sample_merged_bam):
    print("Multiple unmapped bams found; merge requested")
    input_bams = args.unmapped_bam
    print("input bams: " + str(input_bams))
    merge_folder = os.path.join(paths.pipeline_outfolder, "unmapped_bam/")
    input_string = " INPUT=" + " INPUT=".join(input_bams)
    output_merge = os.path.join(merge_folder, sample_merged_bam)
    cmd = "java -jar " + os.path.join(paths.picard_dir, "MergeSamFiles.jar")
    cmd += input_string
    cmd += " OUTPUT=" + output_merge
    cmd += " ASSUME_SORTED=TRUE"
    cmd += " CREATE_INDEX=TRUE"

    pipetk.call_lock(cmd, "lock.merge", paths.pipeline_outfolder, output_merge)
    args.unmapped_bam = sample_merged_bam  #update unmapped bam reference
else:
    # Link the file into the unmapped_bam directory
    print("Single unmapped bam found; no merge required")
    print("Unmapped bam:\t\t" + str(args.unmapped_bam[0]))
    args.unmapped_bam = args.unmapped_bam[0]
    local_unmapped_bam = paths.pipeline_outfolder + "unmapped_bam/" + args.sample_name + ".bam"
    call("ln -s " + args.unmapped_bam + " " + local_unmapped_bam, shell=True)

print("Input Unmapped Bam: " + args.unmapped_bam)
#check for file exists:
if not os.path.isfile(args.unmapped_bam):
    print args.unmapped_bam + "is not a file"

# Fastq conversion
########################################################################################
# Uses ngstk module.
out_fastq = os.path.join(paths.pipeline_outfolder, "fastq/", args.sample_name)
ngstk.bam_to_fastq(args.unmapped_bam, out_fastq, args.paired_end, paths, args.sanity_check)

# Adapter trimming
########################################################################################
pipetk.timestamp("### Adapter trimming: ")

if not args.paired_end:
    cmd = "java -jar  " + paths.trimmomatic_jar + " SE -phred33 -threads 30"
    cmd += " -trimlog " + os.path.join(paths.pipeline_outfolder, "fastq/") + "trimlog.log "
    cmd += out_fastq + "_R1.fastq "
    cmd += out_fastq + "_R1_trimmed.fastq "
    cmd += " HEADCROP:6 ILLUMINACLIP:" + paths.adapter_file + ":2:10:4:1:true:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:25"

else:
    cmd = "java -jar " + paths.trimmomatic_jar + " PE -phred33 -threads 30"
    cmd += " -trimlog " + os.path.join(paths.pipeline_outfolder, "fastq/") + "trimlog.log "
    cmd += out_fastq + "_R1.fastq "
    cmd += out_fastq + "_R2.fastq "
    cmd += out_fastq + "_R1_trimmed.fastq "
    cmd += out_fastq + "_R1_unpaired.fastq "
    cmd += out_fastq + "_R2_trimmed.fastq "
    cmd += out_fastq + "_R2_unpaired.fastq "
    cmd += " HEADCROP:6 ILLUMINACLIP:" + paths.adapter_file + ":2:10:4:1:true:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:25"

pipetk.call_lock(cmd, "lock.adapter", paths.pipeline_outfolder, out_fastq + "_R1_trimmed.fastq")
trimmed_fastq = out_fastq + "_R1_trimmed.fastq"

if (args.sanity_check):
    x = ngstk.count_reads(trimmed_fastq)
    pipetk.report_result("Trimmed size", x, paths)

# RNA BitSeq pipeline.
########################################################################################
pipetk.timestamp("### Bowtie1 alignment: ")
bowtie1_folder = paths.pipeline_outfolder + "/bowtie1_" + args.genome_assembly + "/"
pipetk.make_sure_path_exists(bowtie1_folder)
out_bowtie1 = bowtie1_folder + args.sample_name + ".aln.sam"

if not args.paired_end:
    cmd = paths.bowtie1
    cmd += " -q -p 6 -a -m 100 --sam "
    cmd += paths.bowtie_indexed_genome + " "
    cmd += out_fastq + "_R1_trimmed.fastq"
    cmd += " " + out_bowtie1
else:
    cmd = paths.bowtie1
    cmd += " -q -p 6 -a -m 100 --minins 0 --maxins 5000 --fr --sam --chunkmbs 200 "
    cmd += paths.bowtie_indexed_genome
    cmd += " -1 " + out_fastq + "_R1_trimmed.fastq"
    cmd += " -2 " + out_fastq + "_R2_trimmed.fastq"
    cmd += " " + out_bowtie1

pipetk.call_lock(cmd, "lock.bowtie1", paths.pipeline_outfolder, out_bowtie1)

if (args.sanity_check):
    x = ngstk.count_reads(out_bowtie1 , " -S -F4 " , args.paired_end )
    pipetk.report_result("Aligned reads", x, paths)

pipetk.timestamp("### Raw: SAM to BAM conversion and sorting: ")

cmd = ngstk.sam_conversions(out_bowtie1,False)
pipetk.call_lock(cmd, "lock.sam2bam", paths.pipeline_outfolder, out_bowtie1.replace(".sam" , "_sorted.bam"),shell=True)


pipetk.timestamp("### Aligned read filtering: ")

out_sam_filter = bowtie1_folder + args.sample_name + ".aln.filt.sam"
headerLines = subprocess.check_output("samtools view -SH " + out_bowtie1 + "|wc -l", shell=True).strip()
cmd = "python " + paths.src_folder + "/bisulfiteReadFiltering_forRNA.py"
cmd += " --infile=" + out_bowtie1
cmd += " --outfile=" + out_sam_filter
cmd += " --skipped=" + out_sam_filter.replace(".filt." , ".skipped.")
cmd += " --skipHeaderLines=" + headerLines
cmd += " --genome=" + args.genome_assembly
cmd += " --genomeDir=" + paths.ref_genome
cmd += " --minNonCpgSites=3"
cmd += " --minConversionRate=0.9"
cmd += " --maxConversionRate=0.1"
cmd += " -r"


if args.paired_end:
    cmd = cmd + " --pairedEnd"

pipetk.call_lock(cmd, "lock.read_filter", paths.pipeline_outfolder, out_sam_filter)

if (args.sanity_check):
    x = ngstk.count_reads(out_sam_filter , " -S " , args.paired_end)
    pipetk.report_result("Filtered", x, paths)


pipetk.timestamp("### Set flag in skipped reads to 4 (unmapped): ")
joined_sam = out_sam_filter.replace(".filt." , ".filtFlag.")
skipped_sam = out_sam_filter.replace(".filt." , ".skipped.")
cmd = "samtools view -hS " + out_sam_filter + " > " + joined_sam + "\n"
cmd += "awk -v skip=" + headerLines + " -v OFS=\"\\t\" '{if (NR>skip){$2=4;$3=\"*\";$4=0;$5=0;$6=\"*\";$7=\"*\";$8=0;$9=0; print}}' " + skipped_sam
cmd += " >>" + joined_sam

pipetk.call_lock(cmd, "lock.changeFlag", paths.pipeline_outfolder, joined_sam , shell=True)

if (args.sanity_check):
    x = ngstk.count_reads(joined_sam , " -S " , args.paired_end)
    pipetk.report_result("Filtered + skipped", x, paths)

pipetk.timestamp("### Filtered: SAM to BAM conversion, sorting and depth calculation: ")
cmd = ngstk.sam_conversions(joined_sam)
pipetk.call_lock(cmd, "lock.filt_sam2bam", paths.pipeline_outfolder, joined_sam.replace(".sam" , "_sorted.depth"),shell=True)

pipetk.timestamp("### Skipped: SAM to BAM conversion and sorting: ")
cmd = ngstk.sam_conversions(skipped_sam,False)
pipetk.call_lock(cmd, "lock.skip_sam2bam", paths.pipeline_outfolder, skipped_sam.replace(".sam" , "_sorted.bam"),shell=True)



# BitSeq
########################################################################################
pipetk.timestamp("### Expression analysis (BitSeq): ")

bitSeq_dir = bowtie1_folder + "/bitSeq"
pipetk.make_sure_path_exists(bitSeq_dir)
out_bitSeq = bitSeq_dir + "/" + args.sample_name + ".counts"

cmd = "Rscript " + paths.scripts_dir + "/bitSeq_parallel.R " + " " + joined_sam + " " + bitSeq_dir + " " + paths.ref_genome_fasta

pipetk.call_lock(cmd, "lock.bitSeq", paths.pipeline_outfolder, out_bitSeq)


# ERCC Spike-in alignment
########################################################################################

pipetk.timestamp("### ERCC: Convert unmapped reads into fastq files: ")

unmappable_bam = out_bowtie1.replace(".sam","_unmappable")
cmd = "samtools view -hbS -f4 " + out_bowtie1 + " > " + unmappable_bam + ".bam"
pipetk.call_lock(cmd, "lock.unmappable_bam", paths.pipeline_outfolder, unmappable_bam,shell=True)

ngstk.bam_to_fastq(unmappable_bam + ".bam", unmappable_bam , args.paired_end, paths, args.sanity_check)


pipetk.timestamp("### ERCC: Bowtie1 alignment: ")
bowtie1_folder = paths.pipeline_outfolder + "/bowtie1_" + args.ERCC_assembly + "/"
pipetk.make_sure_path_exists(bowtie1_folder)
out_bowtie1 = bowtie1_folder + args.sample_name + "_ERCC.aln.sam"

if not args.paired_end:
    cmd = paths.bowtie1
    cmd += " -q -p 6 -a -m 100 --sam "
    cmd += paths.bowtie_indexed_ERCC + " "
    cmd += unmappable_bam + "_R1.fastq"
    cmd += " " + out_bowtie1
else:
    cmd = paths.bowtie1
    cmd += " -q -p 6 -a -m 100 --minins 0 --maxins 5000 --fr --sam --chunkmbs 200 "
    cmd += paths.bowtie_indexed_ERCC
    cmd += " -1 " + unmappable_bam + "_R1.fastq"
    cmd += " -2 " + unmappable_bam + "_R2.fastq"
    cmd += " " + out_bowtie1

pipetk.call_lock(cmd, "lock.ERCCbowtie1", paths.pipeline_outfolder, out_bowtie1)

if (args.sanity_check):
    x = ngstk.count_reads(out_bowtie1 , " -S -F 4 " , args.paired_end )
    pipetk.report_result("Aligned reads", x, paths)


pipetk.timestamp("### ERCC: SAM to BAM conversion, sorting and depth calculation: ")
cmd = ngstk.sam_conversions(out_bowtie1)
pipetk.call_lock(cmd, "lock.ERCC_sam2bam", paths.pipeline_outfolder, out_bowtie1.replace(".sam" , "_sorted.depth"),shell=True)


# BitSeq
########################################################################################
pipetk.timestamp("### Expression analysis (BitSeq): ")

bitSeq_dir = bowtie1_folder + "/bitSeq"
pipetk.make_sure_path_exists(bitSeq_dir)
out_bitSeq = bitSeq_dir + "/" + out_bowtie1.replace(".aln.sam" , ".counts")

cmd = "Rscript " + paths.scripts_dir + "/bitSeq_parallel.R " + " " + out_bowtie1 + " " + bitSeq_dir + " " + paths.ref_ERCC_fasta

pipetk.call_lock(cmd, "lock.bitSeq", paths.pipeline_outfolder, out_bitSeq)



# Cleanup
########################################################################################


# remove temporary marker file:
pipetk.stop_pipeline(paths, args, start_time , pipeline )





