
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
		bam_size = subprocess.check_output("samtools view -c " + bam_file, shell=True)
		print ("Bam size: " + str(bam_size))
		fastq_size = subprocess.check_output("wc -l " + out_fastq_pre + "_R1.fastq | cut -f1 -d' '", shell=True)
		if not paired_end:
			fastq_reads = int(fastq_size) / 4
		else:
			fastq_reads = int(fastq_size) / 2

		print ("Fastq reads: " + "{:,}".format(fastq_reads))
		if (fastq_reads!= int(bam_size)):
			raise Exception("Fastq conversion error? Size doesn't match unaligned bam")
