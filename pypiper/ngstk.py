#!/usr/env python

import os
import subprocess
import errno
from AttributeDict import AttributeDict as _AttributeDict

# Path variables
# #######################################################################################
# in here is stuff that will be used by multiple pipelines.


# HOW TO USE IN PIPELINES:
# from pypiper import ngstk
# myngstk = ngstk.NGSTk(config_file)
# myngstk.bam_to_fastq(arguments NOT REQUIRING config.tools OBJECT)

class NGSTk(_AttributeDict):
	"""
	Class to hold functions to build command strings used during pipeline runs.
	Object is instantiated with a string of a path to a yaml `pipeline config file`.
	Since NGSTk inherits from `AttributeDict`, the passed config file and its elements \
	will be accessible through the NGSTk object as attributes under `config` (e.g. \
	`NGSTk.config.tools.java`).

	Example:
	from pypiper.ngstk import NGSTk
	myngstk = NGSTk("pipeline_config_file")
	myngstk.bam_to_fastq(arguments NOT REQUIRING config.tools OBJECT)

	:param config_file: Path to pipeline yaml config file.
	:type config_file: str
	"""

	def __init__(self, config_file):
		# parse yaml into the project's attributes
		# self.add_entries(**config)

		with open(config_file, 'r') as config_file:
			import yaml
			config = yaml.load(config_file)
		super(NGSTk, self).__init__(config, default=True)

	def make_sure_path_exists(self, path):
		try:
			os.makedirs(path)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise

	def markDuplicates(self, aligned_file, out_file, metrics_file, remove_duplicates="True"):
		cmd = self.tools.java + " -Xmx2g -jar " + self.tools.picard + " MarkDuplicates"
		cmd += " INPUT=" + aligned_file
		cmd += " OUTPUT=" + out_file
		cmd += " METRICS_FILE=" + metrics_file
		cmd += " REMOVE_DUPLICATES=" + remove_duplicates
		return cmd

	def bam_to_fastq(self, bam_file, out_fastq_pre, paired_end):
		self.make_sure_path_exists(os.path.dirname(out_fastq_pre))
		cmd = self.tools.java + " -Xmx2g -jar " + self.tools.picard + " SamToFastq"
		cmd += " I=" + bam_file
		cmd += " F=" + out_fastq_pre + "_R1.fastq"
		if paired_end:
			cmd += " F2=" + out_fastq_pre + "_R2.fastq"
		cmd += " INCLUDE_NON_PF_READS=true"
		cmd += " QUIET=true"
		cmd += " VERBOSITY=ERROR"
		return cmd

	def merge_bams(self, input_bams, merged_bam):
		if not len(input_bams) > 1:
			print("No merge required")
			return 0

		print("Merging multiple bams: " + str(input_bams))
		input_string = " INPUT=" + " INPUT=".join(input_bams)
		cmd = self.tools.java + " -Xmx2g -jar " + self.tools.picard + " MergeSamFiles"
		cmd += input_string
		cmd += " OUTPUT=" + merged_bam
		cmd += " ASSUME_SORTED=TRUE"
		cmd += " CREATE_INDEX=TRUE"
		return(cmd)

	def count_lines(self, fileName):
		"""
		Uses the command-line utility wc to count the number of lines in a file.
		:param file: Filename
		"""
		x = subprocess.check_output("wc -l " + fileName + " | cut -f1 -d' '", shell=True)
		return x

	###################################
	# Read counting functions
	###################################
	# In these functions, A paired-end read, with 2 sequences, counts as a two reads

	def count_unique_reads(self, fileName, paired_end):
		"""
		Sometimes alignment software puts multiple locations for a single read; if you just count
		those reads, you will get an inaccurate count. This is _not_ the same as multimapping reads,
		which may or may not be actually duplicated in the bam file (depending on the alignment
		software).
		This function counts each read only once.
		This accounts for paired end or not for free because pairs have the same read name.
		In this function, a paired-end read would count as 2 reads.
		"""
		if fileName.endswith("sam"):
			param = "-S"
		if fileName.endswith("bam"):
			param = ""
		if paired_end:
			r1 = self.samtools_view(fileName, param=param + " -f64", postpend=" | cut -f1 | sort -k1,1 -u | wc -l ")
			r2 = self.samtools_view(fileName, param=param + " -f128", postpend=" | cut -f1 | sort -k1,1 -u | wc -l ")
		else:
			r1 = self.samtools_view(fileName, param=param + "", postpend=" | cut -f1 | sort -k1,1 -u | wc -l ")
			r2 = 0
		return int(r1) + int(r2)

	def count_unique_mapped_reads(self, fileName, paired_end):
		"""
		For a bam or sam file with paired or or single-end reads, returns the
		number of mapped reads, counting each read only once, even if it appears
		mapped at multiple locations.
		:param file: Filename
		:param paired_end: True/False paired end data
		"""
		if fileName.endswith("sam"):
			param = "-S -F4"
		if fileName.endswith("bam"):
			param = "-F4"
		if paired_end:
			r1 = self.samtools_view(fileName, param=param + " -f64", postpend=" | cut -f1 | sort -k1,1 -u | wc -l ")
			r2 = self.samtools_view(fileName, param=param + " -f128", postpend=" | cut -f1 | sort -k1,1 -u | wc -l ")
		else:
			r1 = self.samtools_view(fileName, param=param + "", postpend=" | cut -f1 | sort -k1,1 -u | wc -l ")
			r2 = 0
		return int(r1) + int(r2)

	def count_flag_reads(self, fileName, flag, paired_end):
		"""
		Counts the number of reads with the specified flag.
		:param flag: sam flag value to be read
		:param paired_end: This parameter is ignored; samtools automatically correctly responds depending
		on the data in the bamfile. We leave the option here just for consistency, since all the other
		counting functions require the parameter. This makes it easier to swap counting functions during
		pipeline development.
		"""
		param = " -c -f" + str(flag)
		if fileName.endswith("sam"):
			param += " -S"
		return self.samtools_view(fileName, param=param)

	def count_multimapping_reads(self, fileName, paired_end):
		"""
		Counts the number of reads that mapped to multiple locations. Warning:
		currently, if the alignment software includes the reads at multiple locations, this function
		will count those more than once. This function is for software that randomly assigns,
		but flags reads as multimappers.
		:param paired_end: This parameter is ignored; samtools automatically correctly responds depending
		on the data in the bamfile. We leave the option here just for consistency, since all the other
		counting functions require the parameter. This makes it easier to swap counting functions during
		pipeline development.
		"""
		return int(self.count_flag_reads(fileName, 256, paired_end))

	def count_uniquelymapping_reads(self, fileName, paired_end):
		"""
		Counts the number of reads that mapped to a unique position.
		:param paired_end: This parameter is ignored.
		"""
		param = " -c -F256"
		if fileName.endswith("sam"):
			param += " -S"
		return self.samtools_view(fileName, param=param)

	def count_fail_reads(self, fileName, paired_end):
		"""
		Counts the number of reads that failed platform/vendor quality checks.
		:param paired_end: This parameter is ignored; samtools automatically correctly responds depending
		on the data in the bamfile. We leave the option here just for consistency, since all the other
		counting functions require the parameter. This makes it easier to swap counting functions during
		pipeline development.
		"""
		return int(self.count_flag_reads(fileName, 512, paired_end))

	def samtools_view(self, fileName, param, postpend=""):
		"""
		Runs samtools view on a file, with flexible parameters and post-processing. Used internally to implement
		the various count_reads functions.
		:param file: Filename
		:param param: String of parameters to pass to samtools view
		:param postpend: String to postpend to the samtools command;
		useful to add cut, sort, wc operations to the samtools view output.
		"""
		x = subprocess.check_output(self.tools.samtools + " view " + param + " " + fileName + " " + postpend, shell=True)
		return x

	def count_reads(self, fileName, paired_end):
		"""
		To do: remove the default paired_end option, you should be forced to specify.
		Paired-end reads count as 2 in this function.
		For paired-end reads, this function assumes that the reads are split into 2 fastq files,
		therefore, to count the reads, it divides by 2 instead of 4. This will thus give an incorrect
		result if your paired-end fastq files are in only a single file.
		"""
		if fileName.endswith("bam"):
			return self.samtools_view(fileName, param="-c")
		if fileName.endswith("sam"):
			return self.samtools_view(fileName, param="-c -S")
		if fileName.endswith("fastq") or fileName.endswith("fq"):
			x = self.count_lines(fileName)
			if (paired_end):
				return int(x) / 2
			else:
				return int(x) / 4
		return -1

	def count_mapped_reads(self, fileName, paired_end):
		"""
		Mapped_reads are not in fastq format, so this one doesn't need to accommodate fastq,
		and therefore, doesn't require a paired-end parameter because it only uses samtools view.
		Therefore, it's ok that it has a default parameter, since this is discarded.
		:param pared_end: This parameter is ignored; samtools automatically correctly responds depending
		on the data in the bamfile. We leave the option here just for consistency, since all the other
		counting functions require the parameter. This makes it easier to swap counting functions during
		pipeline development.
		"""
		if fileName.endswith("bam"):
			return self.samtools_view(fileName, param="-c -F4")
		if fileName.endswith("sam"):
			return self.samtools_view(fileName, param="-c -F4 -S")
		return -1

	def sam_conversions(self, sam, depth=True):
		"""
		Convert sam files to bam files, then sort and index them for later use.
		:param depth: also calculate coverage over each position
		"""
		cmd = self.tools.samtools + " view -bS " + sam + " > " + sam.replace(".sam", ".bam") + "\n"
		cmd += self.tools.samtools + " sort " + sam.replace(".sam", ".bam") + " " + sam.replace(".sam", "_sorted") + "\n"
		cmd += self.tools.samtools + " index " + sam.replace(".sam", "_sorted.bam") + "\n"
		if depth:
			cmd += self.tools.samtools + " depth " + sam.replace(".sam", "_sorted.bam") + " > " + sam.replace(".sam", "_sorted.depth") + "\n"
		return cmd

	def bam_conversions(self, bam, depth=True):
		"""
		Sort and index bam files for later use.
		:param depth: also calculate coverage over each position
		"""
		cmd = self.tools.samtools + " view -h " + bam + " > " + bam.replace(".bam", ".sam") + "\n"
		cmd += self.tools.samtools + " sort " + bam + " " + bam.replace(".bam", "_sorted") + "\n"
		cmd += self.tools.samtools + " index " + bam.replace(".bam", "_sorted.bam") + "\n"
		if depth:
			cmd += self.tools.samtools + " depth " + bam.replace(".bam", "_sorted.bam") + " > " + bam.replace(".bam", "_sorted.depth") + "\n"
		return cmd

	def fastqc(self, bam, outputDir):
		"""Call fastqc on a bam file."""
		cmd = self.tools.fastqc + " --noextract --outdir " + outputDir + " " + bam
		return cmd

	def samtools_index(self, bam):
		"""Index a bam file."""
		cmd = self.tools.samtools + " index {0}".format(bam)
		return cmd

	def slurmHeader(self, jobName, output, queue="shortq", ntasks=1, time="10:00:00", cpusPerTask=16, memPerCpu=2000, nodes=1, userMail=""):
		cmd = """    #!/bin/bash
		#SBATCH --partition={0}
		#SBATCH --ntasks={1}
		#SBATCH --time={2}

		#SBATCH --cpus-per-task={3}
		#SBATCH --mem-per-cpu={4}
		#SBATCH --nodes={5}

		#SBATCH --job-name={6}
		#SBATCH --output={7}

		#SBATCH --mail-type=end
		#SBATCH --mail-user={8}

		# Start running the job
		hostname
		date

		""".format(queue, ntasks, time, cpusPerTask, memPerCpu, nodes, jobName, output, userMail)
		return cmd

	def slurmFooter(self):
		cmd = "    date"
		return cmd

	def slurmSubmitJob(self, jobFile):
		import os
		cmd = "sbatch %s" % jobFile
		return os.system(cmd)

	def removeFile(self, fileName):
		cmd = "rm {0}".format(fileName)
		return cmd

	def moveFile(self, old, new):
		cmd = "mv {0} {1}".format(old, new)
		return cmd

	def makeDir(self, directory):
		cmd = "mkdir -p {0}".format(directory)
		return cmd

	def mergeBams(self, inputBams, outputBam):
		cmd = self.tools.java + " -Xmx2g -jar " + self.tools.picard + " MergeSamFiles"
		cmd += " USE_THREADING=TRUE"
		cmd += " " + (" ".join(["INPUT=%s"] * len(inputBams))) % tuple(inputBams)
		cmd += " OUTPUT={0}".format(outputBam)
		return cmd

	def fastQC(self, inputBam, outputDir, sampleName):
		import os
		initial = os.path.splitext(os.path.basename(inputBam))[0]
		cmd1 = self.tools.fastqc + " --noextract --outdir {0} {1}".format(outputDir, inputBam)
		cmd2 = "mv {0}_fastqc.html {1}_fastqc.html".format(os.path.join(outputDir, initial), os.path.join(outputDir, sampleName))
		cmd3 = "mv {0}_fastqc.zip {1}_fastqc.zip".format(os.path.join(outputDir, initial), os.path.join(outputDir, sampleName))
		return [cmd1, cmd2, cmd3]

	def preseq_curve(self, bam_file, output_prefix):
		return """
		preseq c_curve -B -P -o {0}.yield.txt {1}
		""".format(output_prefix, bam_file)

	def preseq_extrapolate(self, bam_file, output_prefix):
		return """
		preseq lc_extrap -v -B -P -e 1e+9 -o {0}.future_yield.txt {1}
		""".format(output_prefix, bam_file)

	def preseq_coverage(self, bam_file, output_prefix):
		return """
		preseq gc_extrap -o {0}.future_coverage.txt {1}
		""".format(output_prefix, bam_file)

	def bam2fastq(self, inputBam, outputFastq, outputFastq2=None, unpairedFastq=None):
		cmd = self.tools.java + " -Xmx4g -jar " + self.tools.picard + " SamToFastq"
		cmd += " INPUT={0}".format(inputBam)
		cmd += " FASTQ={0}".format(outputFastq)
		if outputFastq2 is not None and unpairedFastq is not None:
			cmd += " SECOND_END_FASTQ={0}".format(outputFastq2)
			cmd += " UNPAIRED_FASTQ={0}".format(unpairedFastq)
		return cmd

	def trimmomatic(
		self, inputFastq1, outputFastq1, cpus, adapters, log,
		inputFastq2=None, outputFastq1unpaired=None,
		outputFastq2=None, outputFastq2unpaired=None):

		PE = False if inputFastq2 is None else True
		pe = "PE" if PE else "SE"
		cmd = self.tools.java + " -Xmx4g -jar " + self.tools.trimmomatic
		cmd += " {0} -threads {1} -trimlog {2} {3}".format(pe, cpus, log, inputFastq1)
		if PE:
			cmd += " {0}".format(inputFastq2)
		cmd += " {0}".format(outputFastq1)
		if PE:
			cmd += " {0} {1} {2}".format(outputFastq1unpaired, outputFastq2, outputFastq2unpaired)
		cmd += " ILLUMINACLIP:{0}:1:40:15:8:true".format(adapters)
		cmd += " LEADING:3 TRAILING:3"
		cmd += " SLIDINGWINDOW:4:10"
		cmd += " MINLEN:36"
		return cmd

	def skewer(self, inputFastq1, outputPrefix, outputFastq1, trimLog, cpus, adapters, inputFastq2=None, outputFastq2=None):

		PE = False if inputFastq2 is None else True
		mode = "pe" if PE else "any"
		cmds = list()
		cmd1 = self.tools.java + " --quiet"
		cmd1 += " -f sanger"
		cmd1 += " -t {0}".format(cpus)
		cmd1 += " -m {0}".format(mode)
		cmd1 += " -x {0}".format(adapters)
		cmd1 += " -o {0}".format(outputPrefix)
		cmd1 += " {0}".format(inputFastq1)
		if inputFastq2 is None:
			cmds.append(cmd1)
		else:
			cmd1 += " {0}".format(inputFastq2)
			cmds.append(cmd1)
		if inputFastq2 is None:
			cmd2 = "mv {0} {1}".format(outputPrefix + "-trimmed.fastq", outputFastq1)
			cmds.append(cmd2)
		else:
			cmd2 = "mv {0} {1}".format(outputPrefix + "-trimmed-pair1.fastq", outputFastq1)
			cmds.append(cmd2)
			cmd3 = "mv {0} {1}".format(outputPrefix + "-trimmed-pair2.fastq", outputFastq2)
			cmds.append(cmd3)
		cmd4 = "mv {0} {1}".format(outputPrefix + "-trimmed.log", trimLog)
		cmds.append(cmd4)
		return cmds

	def bowtie2Map(self, inputFastq1, outputBam, log, metrics, genomeIndex, maxInsert, cpus, inputFastq2=None):
		import re
		outputBam = re.sub("\.bam$", "", outputBam)
		# Admits 2000bp-long fragments (--maxins option)
		cmd = self.tools.bowtie2 + " --very-sensitive -p {0}".format(cpus)
		cmd += " -x {0}".format(genomeIndex)
		cmd += " --met-file {0}".format(metrics)
		if inputFastq2 is None:
			cmd += " {0} ".format(inputFastq1)
		else:
			cmd += " --maxins {0}".format(maxInsert)
			cmd += " -1 {0}".format(inputFastq1)
			cmd += " -2 {0}".format(inputFastq2)
		cmd += " 2> {0} | samtools view -S -b - | samtools sort - {1}".format(log, outputBam)
		return cmd

	def topHatMap(self, inputFastq, outDir, genome, transcriptome, cpus):
		# TODO:
		# Allow paired input
		cmd = self.tools.tophat + " --GTF {0} --b2-L 15 --library-type fr-unstranded --mate-inner-dist 120".format(transcriptome)
		cmd += " --max-multihits 100 --no-coverage-search"
		cmd += " --num-threads {0} --output-dir {1} {2} {3}".format(cpus, outDir, genome, inputFastq)
		return cmd

	def picardMarkDuplicates(self, inputBam, outputBam, metricsFile, tempDir="."):
		import re
		transientFile = re.sub("\.bam$", "", outputBam) + ".dups.nosort.bam"
		outputBam = re.sub("\.bam$", "", outputBam)
		cmd1 = self.tools.java + " -Xmx4g -jar  `which MarkDuplicates.jar`"
		cmd1 += " INPUT={0}".format(inputBam)
		cmd1 += " OUTPUT={0}".format(transientFile)
		cmd1 += " METRICS_FILE={0}".format(metricsFile)
		cmd1 += " VALIDATION_STRINGENCY=LENIENT"
		cmd1 += " TMP_DIR={0}".format(tempDir)
		# Sort bam file with marked duplicates
		cmd2 = self.tools.samtools + " sort {0} {1}".format(transientFile, outputBam)
		# Remove transient file
		cmd3 = "if [[ -s {0} ]]; then rm {0}; fi".format(transientFile)
		return [cmd1, cmd2, cmd3]

	def removeDuplicates(self, inputBam, outputBam, cpus=16):
		cmd = self.tools.sambamba + " markdup -t {0} -r {1} {2}".format(cpus, inputBam, outputBam)
		return cmd

	def filterReads(self, inputBam, outputBam, metricsFile, paired=False, cpus=16, Q=30):
		"""
		Remove duplicates, filter for >Q, remove multiple mapping reads.
		For paired-end reads, keep only proper pairs.
		"""
		import re
		nodups = re.sub("\.bam$", "", outputBam) + ".nodups.nofilter.bam"
		cmd1 = self.tools.sambamba + " markdup -t {0} -r --compression-level=0 {1} {2} 2> {3}".format(cpus, inputBam, nodups, metricsFile)
		cmd2 = self.tools.sambamba + ' view -t {0} -f bam --valid'.format(cpus)
		if paired:
			cmd2 += ' -F "not (unmapped or mate_is_unmapped) and proper_pair'
		else:
			cmd2 += ' -F "not unmapped'
		cmd2 += ' and not (secondary_alignment or supplementary) and mapping_quality >= {0}"'.format(Q)
		cmd2 += ' {0} |'.format(nodups)
		cmd2 += self.tools.sambamba + " sort -t {0} /dev/stdin -o {1}".format(cpus, outputBam)
		cmd3 = "if [[ -s {0} ]]; then rm {0}; fi".format(nodups)
		cmd4 = "if [[ -s {0} ]]; then rm {0}; fi".format(nodups + ".bai")
		return [cmd1, cmd2, cmd3, cmd4]

	def shiftReads(self, inputBam, genome, outputBam):
		import re
		outputBam = re.sub("\.bam$", "", outputBam)
		cmd = self.tools.samtools + " view -h {0} |".format(inputBam)
		cmd += " shift_reads.py {0} |".format(genome)
		cmd += " " + self.tools.samtools + " view -S -b - |"
		cmd += " " + self.tools.samtools + " sort - {0}".format(outputBam)
		return cmd

	def sortIndexBam(self, inputBam, outputBam):
		import re
		tmpBam = re.sub("\.bam", ".sorted", inputBam)
		cmd1 = self.tools.samtools + " sort {0} {1}".format(inputBam, tmpBam)
		cmd2 = "mv {0}.bam {1}".format(tmpBam, outputBam)
		cmd3 = self.tools.samtools + " index {0}".format(outputBam)
		return [cmd1, cmd2, cmd3]

	def indexBam(self, inputBam):
		cmd = self.tools.samtools + " index {0}".format(inputBam)
		return cmd

	def peakTools(self, inputBam, output, plot, cpus):
		cmd = self.tools.Rscript + " `which run_spp.R` -rf -savp -savp={0} -s=0:5:500 -c={1} -out={2}".format(plot, inputBam, output)
		return cmd

	def getFragmentSizes(self, bam):
		try:
			import pysam
			import numpy as np
		except:
			return
		fragSizes = list()
		bam = pysam.Samfile(bam, 'rb')
		for read in bam:
			if bam.getrname(read.tid) != "chrM" and read.tlen < 1500:
				fragSizes.append(read.tlen)
		bam.close()
		return np.array(fragSizes)

	def plotInsertSizesFit(self, bam, plot, outputCSV, maxInsert=1500, smallestInsert=30):
		"""
		Heavy inspiration from here:
		https://github.com/dbrg77/ATAC/blob/master/ATAC_seq_read_length_curve_fitting.ipynb
		"""
		try:
			import pysam
			import numpy as np
			import matplotlib.mlab as mlab
			from scipy.optimize import curve_fit
			from scipy.integrate import simps
			import matplotlib.pyplot as plt
		except:
			return
		try:
			import seaborn as sns
			sns.set_style("whitegrid")
		except:
			pass

		def getFragmentSizes(bam, maxInsert=1500):
			fragSizes = list()
			bam = pysam.Samfile(bam, 'rb')
			for i, read in enumerate(bam):
				if read.tlen < maxInsert:
					fragSizes.append(read.tlen)
			bam.close()
			return np.array(fragSizes)

		def mixtureFunction(x, *p):
			"""
			Mixture function to model four gaussian (nucleosomal) and one exponential (nucleosome-free) distributions.
			"""
			m1, s1, w1, m2, s2, w2, m3, s3, w3, m4, s4, w4, q, r = p
			nfr = expo(x, 2.9e-02, 2.8e-02)
			nfr[:smallestInsert] = 0
			return (
				mlab.normpdf(x, m1, s1) * w1 +
				mlab.normpdf(x, m2, s2) * w2 +
				mlab.normpdf(x, m3, s3) * w3 +
				mlab.normpdf(x, m4, s4) * w4 +
				nfr)

		def expo(x, q, r):
			"""
			Exponential function.
			"""
			return q * np.exp(-r * x)

		# get fragment sizes
		fragSizes = getFragmentSizes(bam)
		# bin
		numBins = np.linspace(0, maxInsert, maxInsert + 1)
		y, scatter_x = np.histogram(fragSizes, numBins, density=1)
		# get the mid-point of each bin
		x = (scatter_x[:-1] + scatter_x[1:]) / 2
		# Parameters are empirical, need to check
		paramGuess = [
			200, 50, 0.7,  # gaussians
			400, 50, 0.15,
			600, 50, 0.1,
			800, 55, 0.045,
			2.9e-02, 2.8e-02  # exponential
		]
		popt3, pcov3 = curve_fit(mixtureFunction, x[smallestInsert:], y[smallestInsert:], p0=paramGuess, maxfev=100000)
		m1, s1, w1, m2, s2, w2, m3, s3, w3, m4, s4, w4, q, r = popt3
		# Plot
		plt.figure(figsize=(12, 12))
		# Plot distribution
		plt.hist(fragSizes, numBins, histtype="step", ec="k", normed=1, alpha=0.5)
		# Plot nucleosomal fits
		plt.plot(x, mlab.normpdf(x, m1, s1) * w1, 'r-', lw=1.5, label="1st nucleosome")
		plt.plot(x, mlab.normpdf(x, m2, s2) * w2, 'g-', lw=1.5, label="2nd nucleosome")
		plt.plot(x, mlab.normpdf(x, m3, s3) * w3, 'b-', lw=1.5, label="3rd nucleosome")
		plt.plot(x, mlab.normpdf(x, m4, s4) * w4, 'c-', lw=1.5, label="4th nucleosome")
		# Plot nucleosome-free fit
		nfr = expo(x, 2.9e-02, 2.8e-02)
		nfr[:smallestInsert] = 0
		plt.plot(x, nfr, 'k-', lw=1.5, label="nucleosome-free")
		# Plot sum of fits
		ys = mixtureFunction(x, *popt3)
		plt.plot(x, ys, 'k--', lw=3.5, label="fit sum")
		plt.legend()
		plt.xlabel("Fragment size (bp)")
		plt.ylabel("Density")
		plt.savefig(plot, bbox_inches="tight")
		# Integrate curves and get areas under curve
		areas = [
			["fraction", "area under curve", "max density"],
			["Nucleosome-free fragments", simps(nfr), max(nfr)],
			["1st nucleosome", simps(mlab.normpdf(x, m1, s1) * w1), max(mlab.normpdf(x, m1, s1) * w1)],
			["2nd nucleosome", simps(mlab.normpdf(x, m2, s2) * w1), max(mlab.normpdf(x, m2, s2) * w2)],
			["3rd nucleosome", simps(mlab.normpdf(x, m3, s3) * w1), max(mlab.normpdf(x, m3, s3) * w3)],
			["4th nucleosome", simps(mlab.normpdf(x, m4, s4) * w1), max(mlab.normpdf(x, m4, s4) * w4)]
		]
		try:
			import csv
			with open(outputCSV, "w") as f:
				writer = csv.writer(f)
				writer.writerows(areas)
		except:
			pass

	def bamToBigWig(self, inputBam, outputBigWig, genomeSizes, genome, tagmented=False, normalize=False):
		import os
		import re
		# TODO:
		# addjust fragment length dependent on read size and real fragment size
		# (right now it asssumes 50bp reads with 180bp fragments)
		cmds = list()
		transientFile = os.path.abspath(re.sub("\.bigWig", "", outputBigWig))
		cmd1 = self.tools.bedtools + " bamtobed -i {0} |".format(inputBam)
		if not tagmented:
			cmd1 += " " + self.tools.bedtools + " slop -i stdin -g {0} -s -l 0 -r 130 |".format(genomeSizes)
			cmd1 += " fix_bedfile_genome_boundaries.py {0} |".format(genome)
		cmd1 += " " + self.tools.genomeCoverageBed + " {0}-bg -g {1} -i stdin > {2}.cov".format(
			"-5 " if tagmented else "",
			genomeSizes,
			transientFile
		)
		cmds.append(cmd1)
		if normalize:
			cmds.append("""awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * 1000000; print}}' {0}.cov {0}.cov > {0}.normalized.cov""".format(transientFile))
		cmds.append(self.tools.bedGraphToBigWig + " {0}{1}.cov {2} {3}".format(transientFile, ".normalized" if normalize else "", genomeSizes, outputBigWig))
		# remove tmp files
		cmds.append("if [[ -s {0}.cov ]]; then rm {0}.cov; fi".format(transientFile))
		if normalize:
			cmds.append("if [[ -s {0}.normalized.cov ]]; then rm {0}.normalized.cov; fi".format(transientFile))
		cmds.append("chmod 755 {0}".format(outputBigWig))
		return cmds

	def addTrackToHub(self, sampleName, trackURL, trackHub, colour, fivePrime=""):
		cmd1 = """echo "track type=bigWig name='{0} {1}' description='{0} {1}'""".format(sampleName, fivePrime)
		cmd1 += """ height=32 visibility=full maxHeightPixels=32:32:25 bigDataUrl={0} color={1}" >> {2}""".format(trackURL, colour, trackHub)
		cmd2 = "chmod 755 {0}".format(trackHub)
		return [cmd1, cmd2]

	def linkToTrackHub(self, trackHubURL, fileName, genome):
		import textwrap
		db = "org" if genome == "hg19" else "db"  # different database call for human
		genome = "human" if genome == "hg19" else genome  # change hg19 to human
		html = """
		<html>
			<head>
				<meta http-equiv="refresh" content="0; url=http://genome.ucsc.edu/cgi-bin/hgTracks?"""
		html += """{db}={genome}&hgt.customText={trackHubURL}" />
			</head>
		</html>
		""".format(trackHubURL=trackHubURL, genome=genome, db=db)
		with open(fileName, 'w') as handle:
			handle.write(textwrap.dedent(html))

	def htSeqCount(self, inputBam, gtf, output):
		sam = inputBam.replace("bam", "sam")
		cmd1 = "samtools view {0} > {1}".format(inputBam, sam)
		cmd2 = "htseq-count -f sam -t exon -i transcript_id -m union {0} {1} > {2}".format(sam, gtf, output)
		cmd3 = "rm {0}".format(sam)
		return [cmd1, cmd2, cmd3]

	def kallisto(self, inputFastq, outputDir, outputBam, transcriptomeIndex, cpus, inputFastq2=None, size=180, b=200):
		cmd1 = self.tools.kallisto + " quant --bias --pseudobam -b {0} -l {1} -i {2} -o {3} -t {4}".format(b, size, transcriptomeIndex, outputDir, cpus)
		if inputFastq2 is None:
			cmd1 += " --single {0}".format(inputFastq)
		else:
			cmd1 += " {0} {1}".format(inputFastq, inputFastq2)
		cmd1 += " | " + self.tools.samtools + " view -Sb - > {0}".format(outputBam)
		cmd2 = self.tools.kallisto + " h5dump -o {0} {0}/abundance.h5".format(outputDir)
		return [cmd1, cmd2]

	def genomeWideCoverage(self, inputBam, genomeWindows, output):
		cmd = self.tools.bedtools + " coverage -abam -counts -a {0} -b {1} > {2}".format(inputBam, genomeWindows, output)
		return cmd

	def calculateFRiP(self, inputBam, inputBed, output):
		cmd = "cut -f 1,2,3 {0} |".format(inputBed)
		cmd += self.tools.bedtools + " coverage -counts -abam {0} -b - |".format(inputBam)
		cmd += " awk '{{sum+=$4}} END {{print sum}}' > {0}".format(output)
		return cmd

	def macs2CallPeaks(self, treatmentBam, outputDir, sampleName, genome, controlBam=None, broad=False):
		sizes = {"hg38": 2.7e9, "hg19": 2.7e9, "mm10": 1.87e9, "dr7": 1.412e9}
		if not broad:
			cmd = self.tools.macs2 + " callpeak -t {0}".format(treatmentBam)
			if controlBam is not None:
				cmd += " -c {0}".format(controlBam)
			cmd += " --bw 200 -g {0} -n {0} --outdir {0}".format(sizes[genome], sampleName, outputDir)
			# --fix-bimodal --extsize 180
		else:
			# Parameter setting for broad factors according to Nature Protocols (2012)
			# Vol.7 No.9 1728-1740 doi:10.1038/nprot.2012.101 Protocol (D) for H3K36me3
			cmd = self.tools.macs2 + " callpeak -t {0}".format(treatmentBam)
			if controlBam is not None:
				cmd += " -c {0}".format(controlBam)
			cmd += " --broad --nomodel --extsize 73 --pvalue 1e-3 -g {0} -n {1} --outdir {2}".format(
				sizes[genome], sampleName, outputDir
			)
		return cmd

	def macs2CallPeaksATACSeq(self, treatmentBam, outputDir, sampleName, genome):
		sizes = {"hg38": 2.7e9, "hg19": 2.7e9, "mm10": 1.87e9, "dr7": 1.412e9}
		cmd = self.tools.macs2 + " callpeak -t {0}".format(treatmentBam)
		cmd += " --nomodel --extsize 147 -g {0} -n {1} --outdir {2}".format(sizes[genome], sampleName, outputDir)
		return cmd

	def macs2PlotModel(self, sampleName, outputDir):
		import os
		# run macs r script
		cmd1 = self.tools.Rscript + " {0}/{1}_model.r".format(os.getcwd(), sampleName)
		# move to sample dir
		cmd2 = "mv {0}/{1}_model.pdf {2}/{1}_model.pdf".format(os.getcwd(), sampleName, outputDir)
		return [cmd1, cmd2]

	def sppCallPeaks(self, treatmentBam, controlBam, treatmentName, controlName, outputDir, broad, cpus):
		broad = "TRUE" if broad else "FALSE"
		cmd = self.tools.Rscript + " `which spp_peak_calling.R` {0} {1} {2} {3} {4} {5} {6}".format(
			treatmentBam, controlBam, treatmentName, controlName, broad, cpus, outputDir
		)
		return cmd

	def bamToBed(self, inputBam, outputBed):
		cmd = self.tools.bedtools + " bamtobed -i {0} > {1}".format(inputBam, outputBed)
		return cmd

	def zinbaCallPeaks(self, treatmentBed, controlBed, cpus, tagmented=False):
		fragmentLength = 80 if tagmented else 180
		cmd = self.tools.Rscript + " `which zinba.R` -l {0} -t {1} -c {2}".format(fragmentLength, treatmentBed, controlBed)
		return cmd

	def filterPeaksMappability(self, peaks, alignability, filteredPeaks):
		cmd = self.tools.bedtools + " intersect -wa -u -f 1"
		cmd += " -a {0} -b {1} > {2} ".format(peaks, alignability, filteredPeaks)
		return cmd

	def homerFindMotifs(self, peakFile, genome, outputDir, size=150, length="8,10,12,14,16", n_motifs=12):
		cmd = "findMotifsGenome.pl {0} {1} {2}".format(peakFile, genome, outputDir)
		cmd += " -mask -size {0} -len {1} -S {2}".format(size, length, n_motifs)
		return cmd

	def AnnotatePeaks(self, peakFile, genome, motifFile, outputBed):
		cmd = "annotatePeaks.pl {0} {1} -mask -mscore -m {2} |".format(peakFile, genome, motifFile)
		cmd += "tail -n +2 | cut -f 1,5,22 > {3}".format(outputBed)
		return cmd

	def centerPeaksOnMotifs(self, peakFile, genome, windowWidth, motifFile, outputBed):

		cmd = "annotatePeaks.pl {0} {1} -size {2} -center {3} |".format(peakFile, genome, windowWidth, motifFile)
		cmd += " awk -v OFS='\t' '{print $2, $3, $4, $1, $6, $5}' |"
		cmd += """ awk -v OFS='\t' -F '\t' '{ gsub("0", "+", $6) ; gsub("1", "-", $6) ; print }' |"""
		cmd += " fix_bedfile_genome_boundaries.py {0} | sortBed > {1}".format(genome, outputBed)
		return cmd

	def getReadType(self, bamFile, n=10):
		"""
		Gets the read type (single, paired) and length of bam file.
		:param bamFile: Bam file to determine read attributes.
		:type bamFile: str
		:param n: Number of lines to read from bam file.
		:type n: int
		:returns: tuple of (readType=string, readLength=int).
		:rtype: tuple
		"""
		import subprocess
		from collections import Counter
		try:
			p = subprocess.Popen([self.tools.samtools, 'view', bamFile], stdout=subprocess.PIPE)
			# Count paired alignments
			paired = 0
			readLength = Counter()
			while n > 0:
				line = p.stdout.next().split("\t")
				flag = int(line[1])
				readLength[len(line[9])] += 1
				if 1 & flag:  # check decimal flag contains 1 (paired)
					paired += 1
				n -= 1
			p.kill()
		except:
			print("Cannot read provided bam file.")
			raise
		# Get most abundant read readLength
		readLength = sorted(readLength)[-1]
		# If at least half is paired, return True
		if paired > (n / 2):
			return ("PE", readLength)
		else:
			return ("SE", readLength)

	def parseBowtieStats(self, statsFile):
		"""
		Parses Bowtie2 stats file, returns series with values.
		:param statsFile: Bowtie2 output file with alignment statistics.
		:type statsFile: str
		"""
		import pandas as pd
		import re
		stats = pd.Series(index=["readCount", "unpaired", "unaligned", "unique", "multiple", "alignmentRate"])
		try:
			with open(statsFile) as handle:
				content = handle.readlines()  # list of strings per line
		except:
			return stats
		# total reads
		try:
			line = [i for i in range(len(content)) if " reads; of these:" in content[i]][0]
			stats["readCount"] = re.sub("\D.*", "", content[line])
			if 7 > len(content) > 2:
				line = [i for i in range(len(content)) if "were unpaired; of these:" in content[i]][0]
				stats["unpaired"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
			else:
				line = [i for i in range(len(content)) if "were paired; of these:" in content[i]][0]
				stats["unpaired"] = stats["readCount"] - int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
			line = [i for i in range(len(content)) if "aligned 0 times" in content[i]][0]
			stats["unaligned"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
			line = [i for i in range(len(content)) if "aligned exactly 1 time" in content[i]][0]
			stats["unique"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
			line = [i for i in range(len(content)) if "aligned >1 times" in content[i]][0]
			stats["multiple"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
			line = [i for i in range(len(content)) if "overall alignment rate" in content[i]][0]
			stats["alignmentRate"] = re.sub("\%.*", "", content[line]).strip()
		except IndexError:
			pass
		return stats

	def parseDuplicateStats(self, statsFile):
		"""
		Parses sambamba markdup output, returns series with values.
		:param statsFile: sambamba output file with duplicate statistics.
		:type statsFile: str
		"""
		import pandas as pd
		import re
		series = pd.Series()
		try:
			with open(statsFile) as handle:
				content = handle.readlines()  # list of strings per line
		except:
			return series
		try:
			line = [i for i in range(len(content)) if "single ends (among them " in content[i]][0]
			series["single-ends"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
			line = [i for i in range(len(content)) if " end pairs...   done in " in content[i]][0]
			series["paired-ends"] = re.sub("\D", "", re.sub("\.\.\..*", "", content[line]))
			line = [i for i in range(len(content)) if " duplicates, sorting the list...   done in " in content[i]][0]
			series["duplicates"] = re.sub("\D", "", re.sub("\.\.\..*", "", content[line]))
		except IndexError:
			pass
		return series

	def parseQC(self, sampleName, qcFile):
		"""
		Parses QC table produced by phantompeakqualtools (spp) and returns sample quality metrics.
		:param sampleName: Sample name.
		:type sampleName: str
		:param qcFile: phantompeakqualtools output file sample quality measurements.
		:type qcFile: str
		"""
		import pandas as pd
		series = pd.Series()
		try:
			with open(qcFile) as handle:
				line = handle.readlines()[0].strip().split("\t")  # list of strings per line
			series["NSC"] = line[-3]
			series["RSC"] = line[-2]
			series["qualityTag"] = line[-1]
		except:
			pass
		return series

	def getPeakNumber(self, sample):
		"""
		Counts number of peaks from a sample's peak file.
		:param sample: A Sample object with the "peaks" attribute.
		:type sampleName: pipelines.Sample
		"""
		import subprocess
		import re
		proc = subprocess.Popen(["wc", "-l", sample.peaks], stdout=subprocess.PIPE)
		out, err = proc.communicate()
		sample["peakNumber"] = re.sub("\D.*", "", out)
		return sample

	def getFRiP(self, sample):
		"""
		Calculates the fraction of reads in peaks for a given sample.
		:param sample: A Sample object with the "peaks" attribute.
		:type sampleName: pipelines.Sample
		"""
		import re
		import pandas as pd
		with open(sample.frip, "r") as handle:
			content = handle.readlines()
		readsInPeaks = int(re.sub("\D", "", content[0]))
		mappedReads = sample["readCount"] - sample["unaligned"]
		return pd.Series(readsInPeaks / mappedReads, index="FRiP")
