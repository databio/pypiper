# This script loops through all the samples,
# submits jobs for them.
import csv
import os
import subprocess
import pipetk


def slurm_submit(template, variables_dict, submit_script, submit=False):
	pipetk.make_sure_path_exists(os.path.dirname(submit_script))
	fin = open(template,'r')
	filedata = fin.read()
	fin.close()
	for key,value in variables_dict.items():
		print str(key) + ":" + str(value)
		filedata = filedata.replace(str(key),str(value))

	fout = open(submit_script,'w')
	fout.write(filedata)
	fout.close()

	if submit:
		subprocess.call("sbatch " + submit_script, shell=True)




psa = "../data/projectSampleAnnotation.csv"

# Read in the annotation tables
f = open(psa, 'rb') # opens the csv file



try:
	input_file = csv.DictReader(f)  # creates the reader object
	for row in input_file:   # iterates the rows of the file in orders
		print row["sample_name"]
		# Construct the path to the demultiplexed, unmapped bam input file:
		unmapped_bam = "/fhgfs/groups/lab_bock/jklughammer/projects/otherProjects/CORE-seq/"
		unmapped_bam += row["experimental_category"] + "/CORE/unmapped_bam/"
		unmapped_bam += row["flowcell"] + "_" +row["lane"] + "__" + row["sample_name"] + ".bam"

		if not os.path.isfile(unmapped_bam):
			print "File not found: " + unmapped_bam
		else:
			cmd = "python wgbs_pipeline.py"
			cmd += " -i " + unmapped_bam
			cmd += " -s " + row["sample_name"]
			print cmd

			# Create new dict
			submit_script = "slurm/" + row["sample_name"] + "_wgbs.sub"
			submit_log = "slurm/" + row["sample_name"] + "_wgbs.log"
			slurm_settings = {}

			slurm_settings["VAR_JOBNAME"] = row["sample_name"]
			slurm_settings["VAR_MEM"] = "6000"
			slurm_settings["VAR_CORES"] = "4"
			slurm_settings["VAR_TIME"] = "00:15:00"
			slurm_settings["VAR_PARTITION"] = "develop"
			slurm_settings["VAR_CODE"] = cmd
			slurm_settings["VAR_LOGFILE"] = submit_log

			slurm_submit("slurm_template.sub", slurm_settings, submit_script)


finally:
	f.close()      # closing





# submit merge files job


# submit RNA (convert fastq.)
# submit DNA methyl (convert fastq.)
# Cleanup? (checks for mapped bam, deletes fastq files).


# Submit cluster jobs.


 # here's my sample code:
#python wgbs_pipeline.py -i /fhgfs/groups/lab_bock/jklughammer/projects/otherProjects/CORE-seq/titration/CORE/unmapped_bam/BSF_0131_C5FD6ACXX_8__CORE_K562_500_1_sub.bam -s test --no-checks


#python wgbs_pipeline.py -i /fhgfs/groups/lab_bock/jklughammer/projects/otherProjects/CORE-seq/titration/CORE/unmapped_bam/BSF_0131_C5FD6ACXX_8__CORE_K562_500_1.bam




