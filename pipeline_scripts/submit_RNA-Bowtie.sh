#! /bin/bash
#Author: Johanna Klughammer <jklughammer@cemm.oeaw.ac.at>
#Version: 2014-10-20
#
#Parameters: 
#1. analysis_dir: needs to be absolute path
#2. genome id eg. hg19 needs to be installed under this name in /fhgfs/prod/ngs_resources/genomes/
#
#Assumptions:
#
#--The shell script(s) referenced within this script are assumed to be in the same folder.
scripts=$(dirname $0)

#--Unmapped BAM files are in the following format: <analysis_dir>/unmapped_bam/BSF_0000_FLOWCELLID_lane#sample_name.bam
#  The sample name is inferred from this path by splitting at # and . and taking the second and first element, respectively.




if [ $# -lt 2 ]; then
	echo "Params! 1: analysis_dir - 2: genome [hg19 | mm10 | rn5]"
	exit 1
fi

analysis_dir=$1
genome_id=$2 # e.g. "hg19_cdna"
ERCC_id=$3 # e.g. "ERCC92" or "" if not needed
setup=$4  # ["PE"|"SE"]
readLen=$5  # [eg 51| 100]
ref_genome_fasta=/fhgfs/groups/lab_bock/jklughammer/projects/comparative_epigenomics/genomes/$genome_id/forBowtie1/$genome_id

if [ $ERCC_id = ""]; then
	ERCC_fasta=""
else
	ERCC_fasta=/fhgfs/groups/lab_bock/jklughammer/projects/comparative_epigenomics/genomes/$ERCC_id/forBowtie1/$ERCC_id
fi

if [ $genome_id == "rn5" ]; then
	echo "Ref genome = $genome_id ---> using Ensembl reference genome"
	ref_genome_fasta=/fhgfs/prod/ngs_resources/genomes/$genome_id/ensembl/$genome_id
fi

if [ ! -d $analysis_dir ]; then
	echo "Analysis directory $analysis_dir does not exist!. Exiting!"
	exit 1
fi

logdir=$analysis_dir/log
mkdir -p $logdir

#replace all # by __

for i in $analysis_dir/unmapped_bam/*#*.bam; do mv "$i" "${i/"#"/__}"; done

# Process each unmapped_bam
for unmapped_bam in `ls $analysis_dir/unmapped_bam/*.bam`; do
	#sample=$(basename $unmapped_bam)
	sample=$(basename $unmapped_bam | cut -f1 -d.)
	echo "Sample: $sample"
	if [ -d "$analysis_dir/${genome_id}_$setup/$sample" ]; then
  		echo "Sample $sample has already been analysed."
#>>>>#		continue
	fi
	
	# run on cluster:
	# qsub -N rrbs_$sample -pe onenode 4 -q longtime.q -e $logdir -o $logdir $scripts/rrbs_pipeline.sh $analysis_dir $unmapped_bam $ref_genome_fasta $genome_id

	# run on new cluster (SLURM)

	sbatch --export=NONE --get-user-env=L --job-name=RNA-Bt_$sample --ntasks=1 --cpus-per-task=6 --mem-per-cpu=4000 --partition=shortq --time=12:00:00 -e $logdir/RNA_bt_${setup}_${sample}_%j.err -o $logdir/RNA_bt_${setup}_${sample}_%j.log $scripts/RNA-Bowtie_pipeline.sh $analysis_dir $unmapped_bam $ref_genome_fasta $genome_id $scripts $setup $ERCC_fasta $readLen
	

	# run locally:
	#sh $scripts/rrbs_pipeline.sh $analysis_dir $unmapped_bam $ref_genome_fasta 2>&1 | tee $logdir/$sample.log
done
