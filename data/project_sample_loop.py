# This script loops through all the samples,
# submits jobs for them.

# Read in the annotation tables
f = open(sys.argv[1], 'rb') # opens the csv file
try:
    reader = csv.reader(f)  # creates the reader object
    for row in reader:   # iterates the rows of the file in orders
        print row    # prints each row
finally:
    f.close()      # closing

def somefunction():
    pass




# submit merge files job


# submit RNA (convert fastq.)
# submit DNA methyl (convert fastq.)
# Cleanup? (checks for mapped bam, deletes fastq files).


# Submit cluster jobs.

