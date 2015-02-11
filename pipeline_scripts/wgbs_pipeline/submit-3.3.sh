#! /bin/bash

if [ $# -lt 3 ]; then
	echo "INPUT PARAMS MISSING. 1: Project_name_dir (e.g. BSA_0034_WGBS_singleCell) - 2: analysis_dir(samples group, folder name only, e.g. KBM7-input-optim) - 3: genome [hg19 | mm10] - 4: Sequencing mode [SR | PE]"
	exit 1
fi




# ----- RUN SETTING -----------------------------------------

project_name=$1
samples_group=$2 # 'KBM7-input-optim' 
genome=$3       # "hg19"
pipeline_version='wgbs-bismark-3.3_'$4

base_folder='/fhgfs/groups/lab_bsf/meth-seq'
project_folder=$base_folder/projects/$project_name  # BSA_0033_WGBS_boston    # BSA_0034_WGBS_singleCell
# Sequencing mode: 1=SR, 2=PE
mode=1

if [ $4 == "PE" ]
then 
	mode=2
fi



# logdir=$project_folder/$samples_group/$pipeline_version/$genome/bismark/log
logdir=$project_folder/$samples_group/$pipeline_version/$genome/log
mkdir -p $logdir
mkdir -p $logdir/tmp

# Process each unmapped_bam
for unmapped_bam in `ls $project_folder/$samples_group/$pipeline_version/unmapped_bam/*.all.unmapped.bam`; do
	sample=$(basename $unmapped_bam)
	echo "Sample: $sample"

	# run on cluster:
	#qsub -N wgbs3_$sample -pe onenode 8 -q bsf.q -e $logdir -o $logdir $base_folder/src/wgbs/dev/wgbs-bismark-3.0.sh $unmapped_bam $samples_group $pipeline_version $genome $base_folder $project_folder $mode
	#qsub -N wgbs_$sample -pe onenode 8 -q longtime.q -e $logdir -o $logdir $base_folder/src/wgbs/dev/wgbs-bismark-3.3.sh $unmapped_bam $samples_group $pipeline_version $genome $base_folder $project_folder $mode

	# run on new cluster (SLURM)
	sbatch --export=NONE --get-user-env=L --job-name=wgbs_$sample --ntasks=1 --cpus-per-task=4 --mem-per-cpu=4000 --time=20-23 --partition=longq -e $logdir/err_$sample_%j.log -o $logdir/out_$sample_%j.log $base_folder/src/wgbs/dev/wgbs-bismark-3.3.sh $unmapped_bam $samples_group $pipeline_version $genome $base_folder $project_folder $mode
	
	# run locally:
	#sh $scripts/rrbs_pipeline.sh $analysis_dir $unmapped_bam $ref_genome_fasta 2>&1 | tee $logdir/$sample.log
done
