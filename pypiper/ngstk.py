#!/usr/env python

import pypiper
import os
import subprocess
import errno

# Path variables
# #######################################################################################
# in here is stuff that will be used by multiple pipelines.

def make_sure_path_exists(path):
	try:
		os.makedirs(path)
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise


def markDuplicates(paths, aligned_file, out_file, metrics_file, remove_duplicates="True"):
	cmd = "java -jar "
	cmd += os.path.join(paths.picard_dir, "MarkDuplicates.jar")
	cmd += " INPUT=" + aligned_file
	cmd += " OUTPUT=" + out_file
	cmd += " METRICS_FILE=" + metrics_file
	cmd += " REMOVE_DUPLICATES=" + remove_duplicates
	return cmd


def bam_to_fastq(bam_file, out_fastq_pre, paired_end, paths):
	make_sure_path_exists(os.path.dirname(out_fastq_pre))
	# Build commands:

	cmd = "java -jar "
	cmd += os.path.join(paths.picard_dir, "SamToFastq.jar")
	cmd += " I=" + bam_file
	cmd += " F=" + out_fastq_pre + "_R1.fastq"

	if paired_end:
		cmd += " F2=" + out_fastq_pre + "_R2.fastq"

	cmd += " INCLUDE_NON_PF_READS=true"
	return cmd


def merge_bams(input_bams, merged_bam, paths):
	if not len(input_bams) > 1:
		print "No merge required"
		return 0

	print("Merging multiple bams: " + str(input_bams))
	input_string = " INPUT=" + " INPUT=".join(input_bams)
	cmd = "java -jar " + os.path.join(paths.picard_dir, "MergeSamFiles.jar")
	cmd += input_string
	cmd += " OUTPUT=" + merged_bam
	cmd += " ASSUME_SORTED=TRUE"
	cmd += " CREATE_INDEX=TRUE"
	return(cmd)


def count_lines(file):
	'''
	Uses the command-line utility wc to count the number of lines in a file.
	:param file: Filename
	'''
	x = subprocess.check_output("wc -l " + file + " | cut -f1 -d' '", shell=True)
	return x



	###################################
	# Read counting functions
	###################################
	# In these functions, A paired-end read, with 2 sequences, counts as a two reads

def count_unique_reads(file, paired_end):
	'''
	Sometimes alignment software puts multiple locations for a single read; if you just count
	those reads, you will get an inaccurate count. This is _not_ the same as multimapping reads,
	which may or may not be actually duplicated in the bam file (depending on the alignment
	software).
	This function counts each read only once.
	This accounts for paired end or not for free because pairs have the same read name.
	In this function, a paired-end read would count as 2 reads.
	'''
	if file.endswith("sam"):
		param = "-S"
	if file.endswith("bam"):
		param = ""
	if paired_end:
		r1 = samtools_view(file, param=param + " -f64", postpend=" | cut -f1 | sort -k1,1 -u | wc -l ")
		r2 = samtools_view(file, param=param + " -f128", postpend=" | cut -f1 | sort -k1,1 -u | wc -l ")
	else:
		r1 = samtools_view(file, param=param + "", postpend=" | cut -f1 | sort -k1,1 -u | wc -l ")
		r2 = 0
	return int(r1) + int(r2)


def count_unique_mapped_reads(file, paired_end):
	'''
	For a bam or sam file with paired or or single-end reads, returns the
	number of mapped reads, counting each read only once, even if it appears
	mapped at multiple locations.
	:param file: Filename
	:param paired_end: True/False paired end data
	'''
	if file.endswith("sam"):
		param = "-S -F4"
	if file.endswith("bam"):
		param = "-F4"
	if paired_end:
		r1 = samtools_view(file, param=param + " -f64", postpend=" | cut -f1 | sort -k1,1 -u | wc -l ")
		r2 = samtools_view(file, param=param + " -f128", postpend=" | cut -f1 | sort -k1,1 -u | wc -l ")
	else:
		r1 = samtools_view(file, param=param + "", postpend=" | cut -f1 | sort -k1,1 -u | wc -l ")
		r2 = 0
	return int(r1) + int(r2)

def count_multimapping_reads(file, paired_end=True):
	'''
	Counts the number of reads that mapped to multiple locations. Warning:
	currently, if the alignment software includes the reads at multiple locations, this function
	will count those more than once. This function is for software that randomly assigns,
	but flags reads as multimappers.
	:param pared_end: This parameter is ignored; samtools automatically correctly responds depending
	on the data in the bamfile. We leave the option here just for consistency, since all the other
	counting functions require the parameter. This makes it easier to swap counting functions during
	pipeline development.
	'''
	param = " -c -f256"
	if file.endswith("sam"):
		param += " -S"

	return samtools_view(file, param=param)
	

def samtools_view(file, param, postpend=""):
	'''
	Runs samtools view on a file, with flexible parameters and post-processing. Used internally to implement
	the various count_reads functions.
	:param file: Filename
	:param param: String of parameters to pass to samtools view
	:param postpend: String to postpend to the samtools command; 
	useful to add cut, sort, wc operations to the samtools view output.
	'''
	x = subprocess.check_output("samtools view " + param + " " + file + " " + postpend, shell=True)
	return x


def count_reads(file, paired_end=True):
	'''
	To do: remove the default paired_end option, you should be forced to specify.
	Paired-end reads count as 2 in this function.
	For paired-end reads, this function assumes that the reads are split into 2 fastq files, 
	therefore, to count the reads, it divides by 2 instead of 4. This will thus give an incorrect
	result if your paired-end fastq files are in only a single file.	
	'''
	if file.endswith("bam"):
		return samtools_view(file, param="-c")
	if file.endswith("sam"):
		return samtools_view(file, param="-c -S")
	if file.endswith("fastq") or file.endswith("fq"):
		x = count_lines(file)
		if (paired_end):
			return int(x) / 2
		else:
			return int(x) / 4
	return -1


def count_mapped_reads(file, paired_end=True):
	'''
	Mapped_reads are not in fastq format, so this one doesn't need to accommodate fastq,
	and therefore, doesn't require a paired-end parameter because it only uses samtools view.
	Therefore, it's ok that it has a default parameter, since this is discarded.
	:param pared_end: This parameter is ignored; samtools automatically correctly responds depending
	on the data in the bamfile. We leave the option here just for consistency, since all the other
	counting functions require the parameter. This makes it easier to swap counting functions during
	pipeline development.
	'''
	if file.endswith("bam"):
		return samtools_view(file, param="-c -F4")
	if file.endswith("sam"):
		return samtools_view(file, param="-c -F4 -S")
	return -1


def sam_conversions(sam, depth=True):
	cmd = "samtools view -bS " + sam + " > " + sam.replace(".sam", ".bam") + "\n"
	cmd += "samtools sort " + sam.replace(".sam", ".bam") + " " + sam.replace(".sam", "_sorted") + "\n"
	cmd += "samtools index " + sam.replace(".sam", "_sorted.bam") + "\n"
	if depth:
		cmd += "samtools depth " + sam.replace(".sam", "_sorted.bam") + " > " + sam.replace(".sam",
																							"_sorted.depth") + "\n"
	return cmd


def bam_conversions(bam, depth=True):
	cmd = "samtools view -h " + bam + " > " + bam.replace(".bam", ".sam") + "\n"
	cmd += "samtools sort " + bam + " " + bam.replace(".bam", "_sorted") + "\n"
	cmd += "samtools index " + bam.replace(".bam", "_sorted.bam") + "\n"
	if depth:
		cmd += "samtools depth " + bam.replace(".bam", "_sorted.bam") + " > " + bam.replace(".bam", "_sorted.depth") + "\n"
	return cmd


def fastqc(bam, outputDir):
	"""Call fastqc on a bam file."""
	cmd = "fastqc --noextract --outdir " + outputDir + " " + bam
	return cmd


def samtools_index(bam):
	"""Index a bam file."""
	cmd = "samtools index {0}".format(bam)
	return cmd
