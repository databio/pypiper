#!/usr/env python

import os
import subprocess
import errno
from AttributeDict import AttributeDict as _AttributeDict


class NGSTk(_AttributeDict):
	"""
	Class to hold functions to build command strings used during pipeline runs.
	Object can be instantiated with a string of a path to a yaml `pipeline config file`.
	Since NGSTk inherits from `AttributeDict`, the passed config file and its elements
	will be accessible through the NGSTk object as attributes under `config` (e.g.
	`NGSTk.tools.java`). In case no `config_file` argument is passed, all commands will
	be returned assuming the tool is in the user's $PATH.

	:param config_file: Path to pipeline yaml config file (optional).
	:type config_file: str

	:Example:

		from pypiper.ngstk import NGSTk as tk
		tk = NGSTk()
		tk.samtools_index("sample.bam")
		# returns: samtools index sample.bam

		# Using a configuration file (custom executable location):
		from pypiper.ngstk import NGSTk
		tk = NGSTk("pipeline_config_file.yaml")
		tk.samtools_index("sample.bam")
		# returns: /home/.local/samtools/bin/samtools index sample.bam

	"""

	def __init__(self, config_file=None, pm=None):
		# parse yaml into the project's attributes
		# self.add_entries(**config)

		if config_file is None:
			super(NGSTk, self).__init__({}, default=True)
		else:
			import yaml
			with open(config_file, 'r') as config_file:
				config = yaml.load(config_file)
			super(NGSTk, self).__init__(config, default=True)

		# Keep a link to the pipeline manager, if one is provided.
		# if None is provided, instantiate "tools" and "parameters" with empty AttributeDicts
		# this allows the usage of the same code for a command with and without using a pipeline manager
		if pm is not None:
			self.pm = pm
			if hasattr(pm.config, "tools"):
				self.tools = self.pm.config.tools
			if hasattr(pm.config, "parameters"):
				self.parameters = self.pm.config.parameters
		else:
			self.tools = _AttributeDict(dict(), default=True)
			self.parameters = _AttributeDict(dict(), default=True)

		# If pigz is available, use that. Otherwise, default to gzip.
		if hasattr(self.pm, "cores") and self.pm.cores > 1 and self.check_command("pigz"):
			self.ziptool = "pigz -p " + self.pm.cores
		else:
			self.ziptool = "gzip"

	def make_dir(self, path):
		try:
			os.makedirs(path)
		except OSError as exception:
			if exception.errno != errno.EEXIST:
				raise

	def make_sure_path_exists(self, path):
		""" Alias for make_dir """
		self.make_dir(path)

	# Borrowed from looper
	def check_command(self, command):
		"""
		Check if command can be called.
		"""
		import os

		# Use `command` to see if command is callable, store exit code
		code = os.system("command -v {0} >/dev/null 2>&1 || {{ exit 1; }}".format(command))

		# If exit code is not 0, report which command failed and return False, else return True
		if code != 0:
			print("Command is not callable: {0}".format(command))
			return False
		else:
			return True

	def get_file_size(self, filenames):
		"""
		Get size of all files in string (space-separated) in megabytes (Mb).

		:param filenames: a space-separated string of filenames
		:type filenames: str
		"""
		# use (1024 ** 3) for gigabytes
		# equivalent to: stat -Lc '%s' filename

		# If given a list, recurse through it.
		if type(filenames) is list:
			return(sum([self.get_file_size(filename) for filename in filenames]))

		return round(sum([float(os.stat(f).st_size) for f in filenames.split(" ")]) / (1024 ** 2), 4)

	def markDuplicates(self, aligned_file, out_file, metrics_file, remove_duplicates="True"):
		cmd = self.tools.java
		if self.pm.javamem:  # If a memory restriction exists.
			cmd += " -Xmx" + self.pm.javamem
		cmd += " -jar " + self.tools.picard + " MarkDuplicates"
		cmd += " INPUT=" + aligned_file
		cmd += " OUTPUT=" + out_file
		cmd += " METRICS_FILE=" + metrics_file
		cmd += " REMOVE_DUPLICATES=" + remove_duplicates
		return cmd

	def bam_to_fastq(self, bam_file, out_fastq_pre, paired_end):
		self.make_sure_path_exists(os.path.dirname(out_fastq_pre))
		cmd = self.tools.java + " -Xmx" + self.pm.javamem
		cmd += " -jar " + self.tools.picard + " SamToFastq"
		cmd += " I=" + bam_file
		cmd += " F=" + out_fastq_pre + "_R1.fastq"
		if paired_end:
			cmd += " F2=" + out_fastq_pre + "_R2.fastq"
		cmd += " INCLUDE_NON_PF_READS=true"
		cmd += " QUIET=true"
		cmd += " VERBOSITY=ERROR"
		cmd += " VALIDATION_STRINGENCY=SILENT"
		return cmd

	def bam_to_fastq_awk(self, bam_file, out_fastq_pre, paired_end):
		"""
		This converts bam file to fastq files, but using awk. As of 2016, this is much faster 
		than the standard way of doing this using Picard, and also much faster than the 
		bedtools implementation as well; however, it does no sanity checks and assumes the reads
		are e

		"""
		self.make_sure_path_exists(os.path.dirname(out_fastq_pre))
		fq1 = out_fastq_pre + "_R1.fastq"
		if paired_end:
			fq2 = out_fastq_pre + "_R2.fastq"
			cmd = self.tools.samtools + " view " + bam_file + " | awk '"
			cmd += r'{ if (NR%2==1) print "@"$1"/1\n"$10"\n+\n"$11 > "' + fq1 + '";'
			cmd += r' else print "@"$1"/2\n"$10"\n+\n"$11 > "' + fq2 + '"; }'
			cmd += "'"  # end the awk command
		else:
			fq2 = None
			cmd = self.tools.samtools + " view " + bam_file + " | awk '"
			cmd += r'{ print "@"$1"\n"$10"\n+\n"$11 > "' + fq1 + '"; }'
			cmd += "'"
		return cmd, fq1, fq2

	def bam_to_fastq_bedtools(self, bam_file, out_fastq_pre, paired_end):
		"""
		Converts bam to fastq; A version using bedtools
		"""
		self.make_sure_path_exists(os.path.dirname(out_fastq_pre))
		fq1 = out_fastq_pre + "_R1.fastq"
		fq2 = None
		cmd = self.tools.bedtools + " bamtofastq -i " + bam_file + " -fq " + fq1 + ".fastq"
		if paired_end:
			fq2 = out_fastq_pre + "_R2.fastq"
			cmd += " -fq2 " + fq2

		return cmd, fq1, fq2


	def get_input_ext(self, input_file):
		"""
		Get the extension of the input_file. Assumes you're using either
		.bam or .fastq/.fq or .fastq.gz/.fq.gz.
		"""
		if input_file.endswith(".bam"):
			input_ext = ".bam"
		elif input_file.endswith(".fastq.gz") or input_file.endswith(".fq.gz"):
			input_ext = ".fastq.gz"
		elif input_file.endswith(".fastq") or input_file.endswith(".fq"):
			input_ext = ".fastq"
		else:
			raise NotImplementedError("This pipeline can only deal with .bam, .fastq, or .fastq.gz files")
		return input_ext

	def merge_or_link(self, input_args, raw_folder, local_base="sample"):
		"""
		This function standardizes various input possibilities by converting
		either .bam, .fastq, or .fastq.gz files into a local file; merging those
		if multiple files given.

		:param local_base: Usually the sample name. This (plus file extension) will
		be the name of the local file linked (or merged) by this function.

		:param input_args: This is a list of arguments, each one is a class of
		inputs (which can in turn be a string or a list). Typically, input_args
		is a list with 2 elements: first a list of read1 files; second
		an (optional!) list of read2 files.
		:type input_args: list
		"""
		self.make_sure_path_exists(raw_folder)

		if type(input_args) != list:
			raise Exception("Input must be a list")

		if any(isinstance(i, list) for i in input_args):
			# We have a list of lists. Process each individually.
			local_input_files = list()
			n_input_files = len(filter(bool, input_args))
			print("Number of input file sets:\t\t" + str(n_input_files))

			for input_i, input_arg in enumerate(input_args):
				# Count how many non-null items there are in the list;
				# we only append _R1 (etc.) if there are multiple input files.
				if n_input_files > 1:
					local_base_extended = local_base + "_R" + str(input_i + 1)
				else:
					local_base_extended = local_base
				if input_arg:
					out = self.merge_or_link(input_arg, raw_folder, local_base_extended)

					print("Local input file: " + out)
					# Make sure file exists:
					if not os.path.isfile(out):
						print out + " is not a file"

					local_input_files.append(out)

			return local_input_files

		else:
			# We have a list of individual arguments. Merge them.

			if len(input_args) == 1:
				# Only one argument in this list. A single input file; we just link
				# it, regardless of file type:
				# Pull the value out of the list
				input_arg = input_args[0]
				input_ext = self.get_input_ext(input_arg)

				# Convert to absolute path
				if not os.path.isabs(input_arg):
					input_arg = os.path.abspath(input_arg)

				# Link it to into the raw folder
				local_input_abs = os.path.join(raw_folder, local_base + input_ext)
				self.pm.run(
					"ln -sf " + input_arg + " " + local_input_abs,
					target=local_input_abs,
					shell=True)
				# return the local (linked) filename absolute path
				return local_input_abs

			else:
				# Otherwise, there are multiple inputs.
				# If more than 1 input file is given, then these are to be merged
				# if they are in bam format.
				if all([self.get_input_ext(x) == ".bam" for x in input_args]):
					sample_merged = local_base + ".merged.bam"
					output_merge = os.path.join(raw_folder, sample_merged)
					cmd = self.merge_bams(input_args, output_merge)
					self.pm.run(cmd, output_merge)
					cmd2 = self.validate_bam(output_merge)
					self.pm.run(cmd, output_merge, nofail=True)
					return(output_merge)

				# if multiple fastq
				if all([self.get_input_ext(x) == ".fastq.gz" for x in input_args]):
					sample_merged = local_base + ".merged.fastq"
					sample_merged_gz = local_base + ".merged.fastq.gz"
					output_merge = os.path.join(raw_folder, sample_merged)
					output_merge_gz = os.path.join(raw_folder, sample_merged_gz)
					cmd1 = "zcat " + " ".join(input_args) + " > " + output_merge
					cmd2 = "gzip " + output_merge
					self.pm.run([cmd1, cmd2], output_merge_gz)
					return(output_merge_gz)

				if all([self.get_input_ext(x) == ".fastq" for x in input_args]):
					sample_merged = local_base + ".merged.fastq"
					output_merge = os.path.join(raw_folder, sample_merged)
					cmd = "cat " + " ".join(input_args) + " > " + output_merge
					self.pm.run(cmd, output_merge)
					return(output_merge)

				# At this point, we don't recognize the input file types or they
				# do not match.
				raise NotImplementedError("Input files must be of the same type; and can only merge bam or fastq.")

	def input_to_fastq(
		self, input_file, sample_name,
		paired_end, fastq_folder, output_file=None, multiclass=False):
		"""
		Builds a command to convert input file to fastq, for various inputs.

		Takes either .bam, .fastq.gz, or .fastq input and returns
		commands that will create the .fastq file, regardless of input type.
		This is useful to made your pipeline easily accept any of these input
		types seamlessly, standardizing you to the fastq which is still the
		most common format for adapter trimmers, etc.

		It will place the output fastq file in given `fastq_folder`.

		:param input_file: filename of the input you want to convert to fastq
		:type input_file: string

		:returns: A command (to be run with PipelineManager) that will ensure
		your fastq file exists.
		"""

		fastq_prefix = os.path.join(fastq_folder, sample_name)
		self.make_sure_path_exists(fastq_folder)

		# this expects a list; if it gets a string, convert it to a list.
		if type(input_file) != list:
			input_file = [input_file]

		if len(input_file) > 1:
			cmd = []
			output_file = []
			for in_i, in_arg in enumerate(input_file):
				output = fastq_prefix + "_R" + str(in_i + 1) + ".fastq"
				result_cmd, uf, result_file = (self.input_to_fastq(in_arg, sample_name, paired_end, fastq_folder, output, multiclass=True))
				cmd.append(result_cmd)
				output_file.append(result_file)

		else:
			# There was only 1 input class.
			# Convert back into a string
			input_file = input_file[0]
			if not output_file:
				output_file = fastq_prefix + "_R1.fastq"
			input_ext = self.get_input_ext(input_file)

			if input_ext == ".bam":
				print("Found .bam file")
				cmd = self.bam_to_fastq(input_file, fastq_prefix, paired_end)
				# pm.run(cmd, output_file, follow=check_fastq)
			elif input_ext == ".fastq.gz":
				print("Found .fastq.gz file")
				if paired_end and not multiclass:
					# For paired-end reads in one fastq file, we must split the file into 2.
					cmd = self.tools.python + " -u " + os.path.join(self.tools.scripts_dir, "fastq_split.py")
					cmd += " -i " + input_file
					cmd += " -o " + fastq_prefix
				else:
					# For single-end reads, we just unzip the fastq.gz file.
					# or, paired-end reads that were already split.
					cmd = self.ziptool + " -d -c " + input_file + " > " + output_file
					# a non-shell version
					# cmd1 = "gunzip --force " + input_file
					# cmd2 = "mv " + os.path.splitext(input_file)[0] + " " + output_file
					# cmd = [cmd1, cmd2]
			elif input_ext == ".fastq":
				cmd = "ln -sf " + input_file + " " + output_file
				print("Found .fastq file; no conversion necessary")

		return [cmd, fastq_prefix, output_file]

	def check_fastq(self, input_files, output_files, paired_end):
		"""
		Returns a follow sanity-check function to be run after a fastq conversion.
		Run following a command that will produce the fastq files.

		This function will make sure any input files have the same number of reads as the
		output files.
		"""

		# Define a temporary function which we will return, to be called by the
		# pipeline.
		# Must define default parameters here based on the parameters passed in. This locks
		# these values in place, so that the variables will be defined when this function
		# is called without parameters as a follow function by pm.run.

		# This is AFTER merge, so if there are multiple files it means the
		# files were split into read1/read2; therefore I must divide by number
		# of files for final reads.
		def temp_func(input_files=input_files, output_files=output_files, paired_end=paired_end):
			if type(input_files) != list:
				input_files = [input_files]
			if type(output_files) != list:
				output_files = [output_files]
			print(input_files)
			print(output_files)
			n_input_files = len(filter(bool, input_files))
			raw_reads = sum([int(self.count_reads(input_file, paired_end)) for input_file in input_files]) / n_input_files
			self.pm.report_result("Raw_reads", str(raw_reads))
			fastq_reads = sum([int(self.count_reads(output_file, paired_end)) for output_file in output_files]) / n_input_files

			self.pm.report_result("Fastq_reads", fastq_reads)
			input_ext = self.get_input_ext(input_files[0])
			# We can only assess pass filter reads in bam files with flags.
			if input_ext == "bam":
				fail_filter_reads = self.count_fail_reads(input_file, paired_end)
				pf_reads = int(raw_reads) - int(fail_filter_reads)
				self.pm.report_result("PF_reads", str(pf_reads))
			if fastq_reads != int(raw_reads):
				raise Exception("Fastq conversion error? Number of reads doesn't match unaligned bam")

			return fastq_reads

		return temp_func

	def check_trim(self, trimmed_fastq, trimmed_fastq_R2, paired_end, fastqc_folder=None):
		"""
		Returns a follow function for a trimming step. This will count the trimmed reads,
		and run fastqc on the trimmed fastq files.
		"""

		def temp_func(
			trimmed_fastq=trimmed_fastq,
			trimmed_fastq_R2=trimmed_fastq_R2,
			paired_end=paired_end,
			fastqc_folder=fastqc_folder):
			n_trim = float(self.count_reads(trimmed_fastq, paired_end))
			self.pm.report_result("Trimmed_reads", int(n_trim))
			try:
				rr = float(self.pm.get_stat("Raw_reads"))
				self.pm.report_result("Trim_loss_rate", round((rr - n_trim) * 100 / rr, 2))
			except:
				print("Can't calculate trim loss rate without raw read result.")

			# Also run a fastqc (if installed/requested)
			if fastqc_folder:
				self.make_sure_path_exists(fastqc_folder)
				cmd = self.fastqc(trimmed_fastq, fastqc_folder)
				self.pm.run(cmd, lock_name="trimmed_fastqc", nofail=True)
				if paired_end:
					cmd = self.fastqc(trimmed_fastq_R2, fastqc_folder)
					self.pm.run(cmd, lock_name="trimmed_fastqc_R2", nofail=True)

		return temp_func

	def validate_bam(self, input_bam):
		cmd = self.tools.java + " -Xmx" + self.pm.javamem
		cmd += " -jar " + self.tools.picard + " ValidateSamFile"
		cmd += " INPUT=" + input_bam
		return(cmd)

	def merge_bams(self, input_bams, merged_bam, in_sorted="TRUE", tmp_dir=None):
		"""
		The tmp_dir parameter is important because the default can sometimes fill up on
		poorly configured systems.
		"""
		if not len(input_bams) > 1:
			print("No merge required")
			return 0

		print("Merging multiple bams: " + str(input_bams))
		input_string = " INPUT=" + " INPUT=".join(input_bams)
		cmd = self.tools.java + " -Xmx" + self.pm.javamem
		cmd += " -jar " + self.tools.picard + " MergeSamFiles"
		cmd += input_string
		cmd += " OUTPUT=" + merged_bam
		cmd += " ASSUME_SORTED=" + str(in_sorted)
		cmd += " CREATE_INDEX=TRUE"
		cmd += " VALIDATION_STRINGENCY=SILENT"
		if tmp_dir:
			cmd += " TMP_DIR=" + tmp_dir
		return(cmd)

	def count_lines(self, fileName):
		"""
		Uses the command-line utility wc to count the number of lines in a file.
		:param file: Filename
		"""
		x = subprocess.check_output("wc -l " + fileName + " | cut -f1 -d' '", shell=True)
		return x

	def count_lines_zip(self, fileName):
		"""
		Uses the command-line utility wc to count the number of lines in a file.
		For compressed files.
		:param file: Filename
		"""
		x = subprocess.check_output("zcat " + fileName + " | wc -l | cut -f1 -d' '", shell=True)
		return x

	def get_chrs_from_bam(self, fileName):
		"""
		Uses samtools to grab the chromosomes from the header that are contained
		in this bam file.
		"""
		x = subprocess.check_output(self.tools.samtools + " view -H " + fileName + " | grep '^@SQ' | cut -f2| sed s'/SN://'", shell=True)
		# Chromosomes will be separated by newlines; split into list to return
		return x.split()

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
		Paired-end reads count as 2 in this function.
		For paired-end reads, this function assumes that the reads are split into 2
		fastq files, therefore, to count the reads, it divides by 2 instead of 4.
		This will thus give an incorrect result if your paired-end fastq files
		are in only a single file (you must divide by 2 again).
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
		if fileName.endswith("fastq.gz") or fileName.endswith("fq.gz"):
			x = self.count_lines_zip(fileName)
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
		cmd += self.tools.samtools + " sort " + sam.replace(".sam", ".bam") + " -o " + sam.replace(".sam", "_sorted.bam") + "\n"
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

		cmd += self.tools.samtools + " sort " + bam + " -o " + bam.replace(".bam", "_sorted.bam") + "\n"

		cmd += self.tools.samtools + " index " + bam.replace(".bam", "_sorted.bam") + "\n"
		if depth:
			cmd += self.tools.samtools + " depth " + bam.replace(".bam", "_sorted.bam") + " > " + bam.replace(".bam", "_sorted.depth") + "\n"
		return cmd

	def fastqc(self, file, outdir):
		"""Call fastqc on a bam file (or fastq file, right?).
		# You can find the fastqc help with fastqc --help"""
		cmd = self.tools.fastqc + " --noextract --outdir " + outdir + " " + file
		return cmd

	def fastqc_rename(self, input_bam, output_dir, sample_name):
		import os
		cmds = list()
		initial = os.path.splitext(os.path.basename(input_bam))[0]
		cmd1 = self.tools.fastqc + " --noextract --outdir {0} {1}".format(output_dir, input_bam)
		cmds.append(cmd1)
		if not os.path.exists(os.path.join(output_dir, sample_name + "_fastqc.html")):
			cmd2 = "mv {0}_fastqc.html {1}_fastqc.html".format(os.path.join(output_dir, initial), os.path.join(output_dir, sample_name))
			cmd3 = "mv {0}_fastqc.zip {1}_fastqc.zip".format(os.path.join(output_dir, initial), os.path.join(output_dir, sample_name))
			cmds.append(cmd2)
			cmds.append(cmd3)
		return cmds

	def samtools_index(self, bam):
		"""Index a bam file."""
		cmd = self.tools.samtools + " index {0}".format(bam)
		return cmd

	def slurm_header(
		self, job_name, output, queue="shortq", n_tasks=1, time="10:00:00",
		cpus_per_task=8, mem_per_cpu=2000, nodes=1, user_mail="", mail_type="end"):
		cmd = """		#!/bin/bash
		#SBATCH --partition={0}
		#SBATCH --ntasks={1}
		#SBATCH --time={2}

		#SBATCH --cpus-per-task={3}
		#SBATCH --mem-per-cpu={4}
		#SBATCH --nodes={5}

		#SBATCH --job-name={6}
		#SBATCH --output={7}

		#SBATCH --mail-type={8}
		#SBATCH --mail-user={9}

		# Start running the job
		hostname
		date

		""".format(
			queue, n_tasks, time, cpus_per_task, mem_per_cpu,
			nodes, job_name, output, mail_type, user_mail)

		return cmd

	def slurm_footer(self):
		cmd = "		date"
		return cmd

	def slurm_submit_job(self, jobFile):
		import os
		cmd = "sbatch %s" % jobFile
		return os.system(cmd)

	def remove_file(self, file_name):
		cmd = "rm {0}".format(file_name)
		return cmd

	def move_file(self, old, new):
		cmd = "mv {0} {1}".format(old, new)
		return cmd

	def makeDir(self, directory):
		cmd = "mkdir -p {0}".format(directory)
		return cmd

	# REDUNDANT
	# def merge_bams(self, inputBams, outputBam):
	# 	cmd = self.tools.java + " -Xmx" + self.pm.javamem
	# 	cmd += " -jar " + self.tools.picard + " MergeSamFiles"
	# 	cmd += " USE_THREADING=TRUE"
	# 	cmd += " " + (" ".join(["INPUT=%s"] * len(inputBams))) % tuple(inputBams)
	# 	cmd += " OUTPUT={0}".format(outputBam)
	# 	return cmd

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
		cmd = self.tools.java + " -Xmx" + self.pm.javamem
		cmd += " -jar " + self.tools.picard + " SamToFastq"
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
		cmd = self.tools.java + " -Xmx" + self.pm.javamem
		cmd += " -jar " + self.tools.trimmomatic
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
		cmd1 = self.tools.skewer + " --quiet"
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
		# Admits 2000bp-long fragments (--maxins option)
		cmd = self.tools.bowtie2 + " --very-sensitive --no-discordant -p {0}".format(cpus)
		cmd += " -x {0}".format(genomeIndex)
		cmd += " --met-file {0}".format(metrics)
		if inputFastq2 is None:
			cmd += " {0} ".format(inputFastq1)
		else:
			cmd += " --maxins {0}".format(maxInsert)
			cmd += " -1 {0}".format(inputFastq1)
			cmd += " -2 {0}".format(inputFastq2)
		cmd += " 2> {0} | samtools view -S -b - | samtools sort -o {1} -".format(log, outputBam)
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
		cmd1 = self.tools.java + " -Xmx" + self.pm.javamem
		cmd1 += " -jar  `which MarkDuplicates.jar`"
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

	def get_mitochondrial_reads(self, bam_file, output, cpus=4):
		"""
		"""
		tmp_bam = bam_file + "tmp_rmMe"
		cmd1 = self.tools.sambamba + " index -t {0} {1}".format(cpus, bam_file)
		cmd2 = self.tools.sambamba + " slice {0} chrM | {1} markdup -t 4 /dev/stdin {2} 2> {3}".format(bam_file, self.tools.sambamba, tmp_bam, output)
		cmd3 = "rm {}".format(tmp_bam)
		return [cmd1, cmd2, cmd3]

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
		# import re
		# outputBam = re.sub("\.bam$", "", outputBam)
		cmd = self.tools.samtools + " view -h {0} |".format(inputBam)
		cmd += " shift_reads.py {0} |".format(genome)
		cmd += " " + self.tools.samtools + " view -S -b - |"
		cmd += " " + self.tools.samtools + " sort -o {0} -".format(outputBam)
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
		cmd = self.tools.Rscript + " " + self.tools.spp + " -rf -savp -savp={0} -s=0:5:500 -c={1} -out={2}".format(plot, inputBam, output)
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

	def plot_atacseq_insert_sizes(self, bam, plot, output_csv, max_insert=1500, smallest_insert=30):
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
			import matplotlib
			matplotlib.use('Agg')
			import matplotlib.pyplot as plt
		except:
			print("Necessary Python modules couldn't be loaded.")
			return

		try:
			import seaborn as sns
			sns.set_style("whitegrid")
		except:
			pass

		def getFragmentSizes(bam, max_insert=1500):
			fragSizes = list()

			bam = pysam.Samfile(bam, 'rb')

			for i, read in enumerate(bam):
				if read.tlen < max_insert:
					fragSizes.append(read.tlen)
			bam.close()

			return np.array(fragSizes)

		def mixtureFunction(x, *p):
			"""
			Mixture function to model four gaussian (nucleosomal) and one exponential (nucleosome-free) distributions.
			"""
			m1, s1, w1, m2, s2, w2, m3, s3, w3, m4, s4, w4, q, r = p
			nfr = expo(x, 2.9e-02, 2.8e-02)
			nfr[:smallest_insert] = 0

			return (mlab.normpdf(x, m1, s1) * w1 +
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
		numBins = np.linspace(0, max_insert, max_insert + 1)
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

		try:
			popt3, pcov3 = curve_fit(mixtureFunction, x[smallest_insert:], y[smallest_insert:], p0=paramGuess, maxfev=100000)
		except:
			print("Nucleosomal fit could not be found.")
			return

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
		nfr[:smallest_insert] = 0
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

			with open(output_csv, "w") as f:
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
			cmds.append("""awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * 1000; print}}' {0}.cov {0}.cov | sort -k1,1 -k2,2n > {0}.normalized.cov""".format(transientFile))
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
		cmd = self.tools.bedtools + " coverage -counts -abam {0} -b {1} > {2}".format(inputBam, genomeWindows, output)
		return cmd

	def simple_frip(self, input_bam, input_bed):
		cmd = self.tools.samtools + " view -c -L " + input_bed
		cmd += " " + input_bam
		return cmd

	def calculate_FRiP(self, inputBam, inputBed, output, cpus=4):
		cmd = self.tools.sambamba + " depth region -t {0}".format(cpus)
		cmd += " -L {0}".format(inputBed)
		cmd += " {0}".format(inputBam)
		cmd += " | awk '{{sum+=$5}} END {{print sum}}' > {0}".format(output)
		return cmd

	def macs2CallPeaks(treatmentBams, outputDir, sampleName, genome, controlBams=None, broad=False, paired=False):
		"""
		Use MACS2 to call peaks.
		"""
		sizes = {"hg38": 2.7e9, "hg19": 2.7e9, "mm10": 1.87e9, "dr7": 1.412e9}

		cmd = "macs2 callpeak -t {0}".format(treatmentBams if type(treatmentBams) is str else " ".join(treatmentBams))
		if controlBams is not None:
			cmd += " -c {0}".format(controlBams if type(controlBams) is str else " ".join(controlBams))
		if paired:
			cmd += "-f BAMPE "
		if not broad:
			cmd += " --fix-bimodal --extsize 180 --bw 200"
		else:
			cmd += " --broad --nomodel --extsize 73 --pvalue 1e-3"
		cmd += " -g {0} -n {1} --outdir {2}".format(sizes[genome], sampleName, outputDir)

		return cmd

	def macs2CallPeaksATACSeq(self, treatmentBam, outputDir, sampleName, genome):
		sizes = {"hg38": 2.7e9, "hg19": 2.7e9, "mm10": 1.87e9, "dr7": 1.412e9, "mm9": 2.7e9}
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
