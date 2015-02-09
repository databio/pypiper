#!/usr/bin/python2.7
"""
WGBS pipeline
documentation.
"""

from argparse import ArgumentParser
import os
import csv
import os.path
import platform
import sys
from time import sleep
from datetime import datetime
from subprocess import call


print("Script start time: " + str(datetime.now()));
print("Python version:" + platform.python_version())



# Path variables
########################################################################################
picard_dir = "/fhgfs/groups/lab_bsf/meth-seq/src/tools/picard-tools-1.100"
trimmomatic_jar = "/cm/shared/apps/trimmomatic/0.32/trimmomatic-0.32-epignome.jar"
adapter_file = "$base_folder/adapters/epignome_adapters_2-debug.fa"



parser = ArgumentParser(
	description='WGBS-seq pipeline.'
)

parser.add_argument('-i', '--unmapped-bam', default="/fhgfs/groups/lab_bock/shared/projects/unmapped.bam",
					nargs="+",	dest='unmapped_bam', type=str, help="Unmapped bam file(s))")

parser.add_argument('-s', '--sample-name', default="default",
					dest='sample_name', type=str, help='Sample Name')

parser.add_argument('-r', '--project-root', default="/fhgfs/groups/lab_bock/shared/projects/",
					dest='project_root', type=str,
					help='Directory in which the project will reside. Default=/fhgfs/groups/lab_bock/shared/projects/.')

parser.add_argument('-p', '--paired-end', default="True",
					dest='paired_end', type=str, help='Paired End Mode')



args = parser.parse_args()

print(args.project_root)

merge=False

if (len(args.unmapped_bam) > 1):
	print("Multiple unmapped bams found (comma in input); merge requested")
	input_bams = args.unmapped_bam.split(",")
	print("input bams:"+str(input_bams))
	merge=True
	if (args.sample_name == "default"):
		args.sample_name = "merged";
else:
	print("Single unmapped bam found; no merge required")
	if (args.sample_name == "default"):
		args.sample_name = os.path.splitext(os.path.basename(args.unmapped_bam))[0]

print("Sample name: " + args.sample_name)

# put a marker file:
pipeline_temp_marker = 'WGBS-running.temp'

with open(os.path.join(args.project_root + pipeline_temp_marker), 'w') as fout:
	fout.write('')
fout.close()

# Merging
########################################################################################
# If 2 unmapped bam files are given, then these are to be merged.
# Merge them first.
# Wait on merging
# just make sure the unmapped bam exists; if not, wait for it.

# os.path.splitext

if (merge):
	merge_folder = os.path.join(args.project_root, "unmapped_bam/")

	input_string = " INPUT="+" INPUT=".join(input_bams)

	sample_merged_bam = args.sample_name + ".merged.bam"

	cmd = "java -jar " + os.path.join(picard_dir, "MergeSamFiles.jar") \
		  + input_string \
		  + " OUTPUT=" + os.path.join(merge_folder, sample_merged_bam) \
		  + " ASSUME_SORTED=TRUE" \
		  + " CREATE_INDEX=TRUE"

	print(cmd)
	args.unmapped_bam=sample_merged_bam

# Fastq conversion
########################################################################################
# next, make sure fastq exists;


temp_fastq_file = "temp.fastq"
while os.path.isfile(temp_fastq_file):
	print "Fastq conversion is in process..."
	sleep(60)

print "Fastq conversion has not started. Starting conversion..."
# Do fastq conversion here

if not args.paired_end:
	out_fastq = os.path.join(args.project_root, "fastq/", sample_name+".fastq")
	cmd = "java -jar " \
		+ os.path.join(picard_dir, "SamToFastq.jar") \
		+ " I=" + args.unmapped_bam  \
		+ " F=" + out_fastq \
		+ " INCLUDE_NON_PF_READS=true"

else:
	out_fastq_1 = os.path.join(args.project_root, "fastq/", args.sample_name+"_R1.fastq")
	out_fastq_2 = os.path.join(args.project_root, "fastq/", args.sample_name+"_R2.fastq")
	cmd = "java -jar " \
		+ os.path.join(picard_dir, "SamToFastq.jar") \
		+ " I=" + args.unmapped_bam  \
		+ " F=" + out_fastq_1 \
		+ " F2=" + out_fastq_2 \
		+ " INCLUDE_NON_PF_READS=true"

print cmd
print "Fastq conversion complete. Continuing..."

call(["ls", "-l"])

# Adapter trimming
########################################################################################

if not args.paired_end:
	cmd = "java -jar  " + trimmomatic_jar + " SE -phred33 -threads 30" \
	+ " -trimlog " + os.path.join(args.project_root, "fastq/") + "trimlog.log " \
	+ out_fastq \
	+ out_fastq + "_trimmed.fastq" \
	+ " HEADCROP:6 ILLUMINACLIP:" + adapter_file + ":2:10:4:1:true:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16"




else:
	java -jar  $NGS_TRIMMOMATIC/trimmomatic-0.32-epignome.jar PE -phred33 -threads 30 -trimlog $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimming.log $samples_fastq_dir/$(basename $unmapped_bam .bam)_R1.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_R2.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed_R1.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_unpaired_R1.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed_R2.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_unpaired_R2.fastq HEADCROP:6 ILLUMINACLIP:$adapter_file:2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16




# WGBS pipeline.
########################################################################################


# remove temporary marker file:
os.remove(os.path.join(args.project_root + pipeline_temp_marker))

# else:
# cleanup();



print("Script end time: " + str(datetime.now()));






