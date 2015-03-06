import pipetk
import os
import subprocess
# Path variables
# #######################################################################################
# in here is stuff that will be used by multiple pipelines.

def bam_to_fastq(bam_file, out_fastq_pre, paired_end, paths, sanity_check=True):
	pipetk.timestamp("### Fastq conversion: ")
	pipetk.make_sure_path_exists(os.path.dirname(out_fastq_pre))
	# Build commands:

	cmd = "java -jar "
	cmd += os.path.join(paths.picard_dir, "SamToFastq.jar")
	cmd += " I=" + bam_file
	cmd += " F=" + out_fastq_pre + "_R1.fastq"

	if paired_end:
		cmd += " F2=" + out_fastq_pre + "_R2.fastq"

	cmd += " INCLUDE_NON_PF_READS=true"

	pipetk.call_lock(cmd, "lock.fastq", paths.pipeline_outfolder, out_fastq_pre + "_R1.fastq")

	# Sanity checks:
	if (sanity_check):
		bam_size = count_reads(bam_file)
		pipetk.report_result("Bam reads", str(bam_size), paths)
		fastq_size = subprocess.check_output("wc -l " + out_fastq_pre + "_R1.fastq | cut -f1 -d' '", shell=True)
		if not paired_end:
			fastq_reads = int(fastq_size) / 4
		else:
			fastq_reads = int(fastq_size) / 2
		pipetk.report_result("Fastq reads", fastq_reads, paths)
		if (fastq_reads != int(bam_size)):
			raise Exception("Fastq conversion error? Size doesn't match unaligned bam")


def merge_bams(unmapped_bams, merged_bam, paths, sanity_check=True):
	if (len(args.unmapped_bam) > 1):
		merge = True
		if (args.sample_name == "default"):
			args.sample_name = "merged";
	else:
		if (args.sample_name == "default"):
			args.sample_name = os.path.splitext(os.path.basename(args.unmapped_bam[0]))[0]

	if merge and not os.path.isfile(sample_merged_bam):
		print("Multiple unmapped bams found; merge requested")
		input_bams = args.unmapped_bam;
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
		args.unmapped_bam = paths.pipeline_outfolder+"unmapped_bam/" + sample_merged_bam  #update unmapped bam reference
		local_unmapped_bam = paths.pipeline_outfolder+"unmapped_bam/" + sample_merged_bam
	else:
		# Link the file into the unmapped_bam directory
		print("Single unmapped bam found; no merge required")
		print("Unmapped bam:\t\t" + str(args.unmapped_bam[0]))
		args.unmapped_bam = args.unmapped_bam[0]
		local_unmapped_bam = paths.pipeline_outfolder+"unmapped_bam/"+args.sample_name+".bam"
		call("ln -s " + args.unmapped_bam + " " + local_unmapped_bam, shell=True)




def count_lines(file):
	x = subprocess.check_output("wc -l " + file + " | cut -f1 -d' '", shell=True)
	return x


def samtools_view(file, param):
	x = subprocess.check_output("samtools view " + param + " " + file, shell=True)
	return x


def count_reads(file, paired_end=True):
	if file.endswith("bam"):
		return samtools_view(file, param="-c")
	if file.endswith("fastq") or file.endswith("fq"):
		x = count_lines(file)
		if (paired_end):
			return int(x) / 2
		else:
			return int(x) / 4
	if file.endswith("sam"):
		return samtools_view(file, param="-c -S")
	return -1


def count_mapped_reads(file):
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
