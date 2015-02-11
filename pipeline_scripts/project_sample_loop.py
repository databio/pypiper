# This script loops through all the samples,
# submits jobs for them.
import csv
psa = "projectSampleAnnotation.csv"

import pandas as pd
# Read in the annotation tables
f = open(psa, 'rb') # opens the csv file
try:
    reader = csv.reader(f)  # creates the reader object
    for row in reader:   # iterates the rows of the file in orders
        print row    # prints each row
finally:
    f.close()      # closing

def somefunction():
    pass



wgbs_pipeline.py -i unmapped_bam -r COREseq/data/samples --paired-end -g hg19
# submit merge files job


# submit RNA (convert fastq.)
# submit DNA methyl (convert fastq.)
# Cleanup? (checks for mapped bam, deletes fastq files).


# Submit cluster jobs.


 # here's my sample code:

python wgbs_pipeline.py -i /fhgfs/groups/lab_bock/jklughammer/projects/otherProjects/CORE-seq/titration/CORE/unmapped_bam/BSF_0131_C5FD6ACXX_8__CORE_K562_500_1_sub.bam -s test --no-checks












