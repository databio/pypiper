#!/usr/bin/python2.7
"""
WGBS pipeline
documentation.
"""

from argparse import ArgumentParser
import pipetk
import os
import csv
import os.path
import sys
from subprocess import call
from subprocess import check_output
import subprocess

import ngstk

from datetime import datetime

# Argument Parsing
########################################################################################
parser = ArgumentParser(description='WGBS-seq pipeline.')

parser.add_argument('-i', '--unmapped-bam',
					default="/fhgfs/groups/lab_bock/jklughammer/projects/otherProjects/CORE-seq/titration/CORE/unmapped_bam/BSF_0131_C5FD6ACXX_8__CORE_K562_500_1_sub.bam",
					nargs="+", dest='unmapped_bam', help="Input unmapped bam file(s))")

parser.add_argument('-s', '--sample-name', default="default",
					dest='sample_name', type=str, help='Sample Name')

parser.add_argument('-r', '--project-root', default="/fhgfs/groups/lab_bock/nsheffield/COREseq/data/samples/",
					dest='project_root', type=str,
					help='Directory in which the project will reside. Default=/fhgfs/groups/lab_bock/nsheffield/share/COREseq/')

parser.add_argument('-g', '--genome', default="hg19",
					dest='genome_assembly', type=str, help='Genome Assembly')

parser.add_argument('--no-checks', dest='sanity_check', action='store_false', default=True)
parser.add_argument('-p', '--paired-end', dest='paired_end', action='store_true', default=True, help='Paired End Mode')
parser.add_argument('-q', '--single-end', dest='paired_end', action='store_false', help='Single End Mode')

args = parser.parse_args()


# Set up environment path variables
########################################################################################
# Set up an container class to hold paths
class Container:
	pass

paths = Container()
paths.picard_dir = "/fhgfs/groups/lab_bsf/meth-seq/src/tools/picard-tools-1.100"
paths.trimmomatic_jar = "/cm/shared/apps/trimmomatic/0.32/trimmomatic-0.32-epignome.jar"
paths.adapter_file = args.project_root + "data/adapters/epignome_adapters_2_add.fa"
paths.base_folder = "/fhgfs/groups/lab_bsf/meth-seq/"
paths.bismark_indexed_genome = paths.base_folder + "/bisulfite-genomes/bowtie2-indexed/" + args.genome_assembly
paths.src_folder = paths.base_folder + "/src/wgbs/production"
paths.ref_genome = paths.base_folder + "/genomes"
paths.pipeline_outfolder = os.path.join(args.project_root + args.sample_name + "/")
paths.log_file = paths.pipeline_outfolder + "wgbs.log.md"


# Run some initial setting and logging code to start the pipeline.
start_time = pipetk.start_pipeline(paths, args)

# Merging
########################################################################################
# If 2 unmapped bam files are given, then these are to be merged.

merge = False
#if (isinstance(args.unmapped_bam, list)):
if (len(args.unmapped_bam) > 1):
	print("Multiple unmapped bams found; merge requested")
	input_bams = args.unmapped_bam;
	print("input bams: " + str(input_bams))
	merge = True
	if (args.sample_name == "default"):
		args.sample_name = "merged";
else:
	print("Single unmapped bam found; no merge required")
	if (args.sample_name == "default"):
		args.sample_name = os.path.splitext(os.path.basename(args.unmapped_bam))[0]

print("N input bams:\t\t" + str(len(args.unmapped_bam)))
print("Sample name:\t\t" + args.sample_name)


sample_merged_bam = args.sample_name + ".merged.bam"
pipetk.make_sure_path_exists(paths.pipeline_outfolder + "unmapped_bam/")

if merge and not os.path.isfile(sample_merged_bam):
	merge_folder = os.path.join(paths.pipeline_outfolder, "unmapped_bam/")
	input_string = " INPUT=" + " INPUT=".join(input_bams)
	output_merge = os.path.join(merge_folder, sample_merged_bam)
	cmd = "java -jar " + os.path.join(picard_dir, "MergeSamFiles.jar")
	cmd += input_string
	cmd += " OUTPUT=" + output_merge
	cmd += " ASSUME_SORTED=TRUE"
	cmd += " CREATE_INDEX=TRUE"

	pipetk.call_lock(cmd, "lock.merge", paths.pipeline_outfolder, output_merge)
	args.unmapped_bam = sample_merged_bam  #update unmapped bam reference
else:
	# Link the file into the unmapped_bam directory
	args.unmapped_bam = args.unmapped_bam[0]
	local_unmapped_bam = paths.pipeline_outfolder+"unmapped_bam/"+args.sample_name+".bam"
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
	cmd += " HEADCROP:6 ILLUMINACLIP:" + paths.adapter_file + ":2:10:4:1:true:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16"

else:
	cmd = "java -jar " + paths.trimmomatic_jar + " PE -phred33 -threads 30"
	cmd += " -trimlog " + os.path.join(paths.pipeline_outfolder, "fastq/") + "trimlog.log "
	cmd += out_fastq + "_R1.fastq "
	cmd += out_fastq + "_R2.fastq "
	cmd += out_fastq + "_R1_trimmed.fastq "
	cmd += out_fastq + "_R1_unpaired.fastq "
	cmd += out_fastq + "_R2_trimmed.fastq "
	cmd += out_fastq + "_R2_unpaired.fastq "
	cmd += " HEADCROP:6 ILLUMINACLIP:" + paths.adapter_file + ":2:10:4:1:true:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16"

pipetk.call_lock(cmd, "lock.adapter", paths.pipeline_outfolder, out_fastq + "_R1_trimmed.fastq")
trimmed_fastq = out_fastq + "_R1_trimmed.fastq"

# WGBS pipeline.
########################################################################################
pipetk.timestamp("### Bismark alignment: ")
bismark_folder = paths.pipeline_outfolder + "/bismark_" + args.genome_assembly + "/"
pipetk.make_sure_path_exists(bismark_folder)
bismark_temp = bismark_folder + "/" + "bismark_temp"
pipetk.make_sure_path_exists(bismark_temp)
out_bismark = bismark_folder + args.sample_name + ".aln.bam"

if not args.paired_end:
	out_bismark_temp = bismark_folder + args.sample_name + "_R1_trimmed.fastq_bismark_bt2.bam"
	cmd = "bismark " + paths.bismark_indexed_genome + " "
	cmd += out_fastq + "_R1_trimmed.fastq"
	cmd += " --bam --unmapped"
	cmd += " --path_to_bowtie /cm/shared/apps/bowtie/2.2.3/bin --bowtie2 "
	cmd += "--temp_dir " + bismark_temp
	cmd += "--output_dir " + bismark_folder
	cmd += " -p 4 "
	cmd += "; mv " + out_bismark_temp + " " + out_bismark  # mv annoying bismark output
else:
	out_bismark_temp = bismark_folder + args.sample_name + "_R1_trimmed.fastq_bismark_bt2_pe.bam"
	cmd = "bismark " + paths.bismark_indexed_genome
	cmd += " --1 " + out_fastq + "_R1_trimmed.fastq"
	cmd += " --2 " + out_fastq + "_R2_trimmed.fastq"
	cmd += " --bam --unmapped"
	cmd += " --path_to_bowtie /cm/shared/apps/bowtie/2.2.3/bin --bowtie2"
	cmd += " --temp_dir " + bismark_temp
	cmd += " --output_dir " + bismark_folder
	cmd += " --minins 0 --maxins 5000"
	cmd += " -p 4 "
	cmd += "; mv " + out_bismark_temp + " " + out_bismark  # mv annoying bismark output

pipetk.call_lock(cmd, "lock.bismark", paths.pipeline_outfolder, out_bismark)

pipetk.timestamp("### PCR duplicate removal: ")
# Bismark's deduplication forces output naming, how annoying.
out_dedup = bismark_folder + args.sample_name + ".aln.deduplicated.bam"

if not args.paired_end:
	cmd = "deduplicate_bismark --single "
	cmd += out_bismark
	cmd += " --bam"
else:
	cmd = "deduplicate_bismark --paired "
	cmd += out_bismark
	cmd += " --bam"

pipetk.call_lock(cmd, "lock.deduplicate", paths.pipeline_outfolder, out_dedup)

pipetk.timestamp("### Aligned read filtering: ")

# convert bam file into sam file and sort again to
# compensate for a sorting issue of "deduplicate_bismark"
sam_temp = bismark_folder + "/" + "sam_temp"
pipetk.make_sure_path_exists(sam_temp)
out_sam = bismark_folder + args.sample_name + ".aln.deduplicated.sam"
cmd = "samtools view " + out_dedup + " | LANG=C sort "
cmd += " --temporary-directory=" + sam_temp
cmd += " > " + out_sam

pipetk.call_lock(cmd, "lock.sam", paths.pipeline_outfolder, out_sam)

out_sam_filter = bismark_folder + args.sample_name + ".aln.dedup.filt.sam"
cmd = "python " + paths.src_folder + "/bisulfiteReadFiltering.py"
cmd += " --infile=" + out_sam
cmd += " --outfile=" + out_sam_filter
cmd += " --genome=" + args.genome_assembly
cmd += " --genomeDir=" + paths.ref_genome
cmd += " --minNonCpgSites=3"
cmd += " --minConversionRate=0.9"

if args.paired_end:
	cmd = cmd + " --pairedEnd"

pipetk.call_lock(cmd, "lock.read_filter", paths.pipeline_outfolder, out_sam_filter)

pipetk.timestamp("### Methylation calling (bismark extractor): ")

extract_dir = bismark_folder + "extractor"
pipetk.make_sure_path_exists(extract_dir)
out_extractor = extract_dir + "/meth-output.txt"
if not args.paired_end:
	cmd = "bismark_methylation_extractor --single-end"
	cmd += " --report"
	cmd += " --bedGraph"
	cmd += " --merge_non_CpG"
	#cmd += " --cytosine_report"
	cmd += " --genome_folder " + paths.bismark_indexed_genome
	cmd += " --gzip"
	cmd += " --output " + extract_dir
	cmd += " " + out_sam_filter #input file
	cmd += " > " + extract_dir + "/meth-output.txt"
else:
	cmd = "bismark_methylation_extractor --paired-end  --no_overlap"
	cmd += " --report"
	cmd += " --bedGraph"
	cmd += " --merge_non_CpG"
	#cmd += " --cytosine_report"
	cmd += " --genome_folder " + paths.bismark_indexed_genome
	cmd += " --gzip"
	cmd += " --output " + extract_dir
	cmd += " " + out_sam_filter #input file
	cmd += " > " + extract_dir + "/meth-output.txt"

pipetk.call_lock(cmd, "lock.extractor", paths.pipeline_outfolder, out_extractor)


# Final sorting and indexing
########################################################################################
# create sorted and indexed BAM files for visualization and analysis
pipetk.timestamp("### Final sorting and indexing: ")


out_header = bismark_folder + args.sample_name + ".reheader.bam"
out_final = bismark_folder + args.sample_name + ".final.bam"

cmd = "java -jar " + paths.picard_dir + "/ReplaceSamHeader.jar"
cmd += " I=" + out_sam_filter
cmd += " HEADER=" + out_dedup
cmd += " O=" + out_header

pipetk.call_lock(cmd, "lock.ReplaceSamHeader", paths.pipeline_outfolder, output_file=None)


# Sort
cmd = "java -jar " + paths.picard_dir + "/SortSam.jar"
cmd +=" I=" + out_header
cmd +=" O=" + out_final
cmd +=" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
cmd +=" CREATE_INDEX=true"
pipetk.call_lock(cmd, "lock.SortSam", paths.pipeline_outfolder, output_file=None)


# Build index
#cmd = "java -jar " + picard_dir + "/BuildBamIndex.jar"
#cmd += "I=" + out_final
#pipetk.call_lock(cmd, "lock.BuildBamIndex", output_file=None, folder=pipeline_outfolder)

# Cleanup
########################################################################################

os.remove(out_header)


# remove temporary marker file:
pipetk.stop_pipeline(paths, args, start_time)

# remove temporary folders
os.rmdir(bismark_temp)
os.rmdir(sam_temp)


