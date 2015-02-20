
import pipetk
import os
import subprocess
# Path variables
########################################################################################
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
		pipetk.report_result("Fastq reads" , fastq_reads, paths)
		if (fastq_reads!= int(bam_size)):
			raise Exception("Fastq conversion error? Size doesn't match unaligned bam")



def count_lines(file):
	x = subprocess.check_output("wc -l " + file + " | cut -f1 -d' '", shell=True)
	return x


def count_reads_bam(file , param):
	x = subprocess.check_output("samtools view -c " + param + file, shell=True)
	return x


def count_reads(file, param="" , paired_end=True):
	if file.endswith("bam"):
		return count_reads_bam(file , param)
	if file.endswith("fastq") or file.endswith("fq"):
		x = count_lines(file)
		if (paired_end):
			return int(x)/2
		else:
			return int(x)/4
	if file.endswith("sam"):
			return count_reads_bam(file , param)
	return -1


def sam_conversions(sam, depth=True):
    cmd = "samtools view -bS " + sam + " > " + sam.replace(".sam",".bam") + "\n"
    cmd += "samtools sort " + sam.replace(".sam",".bam") + " " + sam.replace(".sam" , "_sorted") + "\n"
    cmd += "samtools index " + sam.replace(".sam" , "_sorted.bam") + "\n"
    if depth:
        cmd += "samtools depth " + sam.replace(".sam" , "_sorted.bam") + " > " + sam.replace(".sam" , "_sorted.depth") + "\n"
    return cmd
