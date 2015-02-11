#! /bin/bash

#Author: Angelo Nuzzo <anuzzo@cemm.oeaw.ac.at>
#Version: 2014-02-03

# ------------------------------------------------------
#  RELEASE NOTE
#
# This version (3.3) is an adaptation of 3.0 including new TRIMMOMATIC version, manually adapted by Johanna & Andreas.
# To be submitted with companion submit-3.3.sh file.

export PYTHONPATH=/cm/shared/apps/cutadapt/1.4.2/lib/python2.7/site-packages:/cm/shared/apps/RSeQC/2.3.9/lib/python2.7/site-packages:/fhgfs/groups/lab_bsf/rrbs_optim/lib/python
export NGS_PICARD=/fhgfs/groups/lab_bsf/meth-seq/src/tools/picard-tools-1.100


if [ $# -eq 0 ]; then
	echo "Params missing! 1: sample unmapped bam file"
	exit 1
fi

unmapped_bam=$1



# ----- RUN SETTING FROM SUBMIT SCRIPT-----------------------------------------

samples_group=$2
pipeline_version=$3
genome=$4
base_folder=$5
project_folder=$6
mode=$7


# ---- PIPELINE-SPECIFIC DIR AND FILES SETTINGS --------------------------------


src_folder=$base_folder/src/wgbs/production

# dir linking to unmapped bam files 
samples_bam_dir=$project_folder/$samples_group/$pipeline_version/unmapped_bam
#mkdir -p $samples_bam_dir

# dir which will contain unampped fastq files
samples_fastq_dir=$project_folder/$samples_group/$pipeline_version/fastq
mkdir -p $samples_fastq_dir
chmod 775 $samples_fastq_dir

# output dir for bismark alignemnt
samples_aligned_dir=$project_folder/$samples_group/$pipeline_version/$genome/bismark
mkdir -p $samples_aligned_dir
chmod 775 $samples_aligned_dir

# output dir for bismark metilation extractor 
meth_dir=$samples_aligned_dir/extractor
mkdir -p $meth_dir
chmod 775 $meth_dir


# Reference genome, bisulfite-indexed by bismark
#bismark_indexed_genome=$project_folder/bisulfite-genome/$genome/   #--Bowtie 1--
bismark_indexed_genome=$base_folder/bisulfite-genomes/bowtie2-indexed/$genome

# Reference (non-converted) genome for python scripts
ref_genome=$base_folder/genomes

# Trimmomatic adapter file
adapter_file=$base_folder/adapters/epignome_adapters_2-debug.fa


 export PYTHONPATH=~/lib/python/



# --------- RENAME UNMAPPED BAM FILES --------- 

#rename .bam .all.unmapped.bam $samples_bam_dir/*.bam 



# -----------------------------------------------
# 1.	CREATE FASTQ + ADAPTER TRIMMING
# -----------------------------------------------

# create fastq files from orginial bam files

	# SINGLE END MODE
	if [ $mode == "1" ]
	then
 	# quick check to skip this step if already done	
		if [ ! -f $samples_fastq_dir/$(basename $unmapped_bam .bam).fastq ]
		then
			echo "Start converting $unmapped_bam to $samples_fastq_dir/$(basename $unmapped_bam .bam).fastq"
			#convert bam to fastq	
			java -jar $NGS_PICARD/SamToFastq.jar I=$unmapped_bam F=$samples_fastq_dir/$(basename $unmapped_bam .bam).fastq INCLUDE_NON_PF_READS=true
	
			# Trim first 6 base from each read, using Trimmomatic (HEADCROP:6)
			java -jar  $NGS_TRIMMOMATIC/trimmomatic-0.32-epignome.jar SE -phred33 -threads 30 -trimlog $samples_fastq_dir/r.log $samples_fastq_dir/$(basename $unmapped_bam .bam).fastq  $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed.fastq  HEADCROP:6 ILLUMINACLIP:$adapter_file:2:10:4:1:true:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16 
			#java -jar  $NGS_TRIMMOMATIC/trimmomatic-0.32.jar SE -phred33 -threads 30 -trimlog $samples_fastq_dir/r.log $samples_fastq_dir/$(basename $unmapped_bam .bam).fastq  $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed.fastq  HEADCROP:6 ILLUMINACLIP:$adapter_file:2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16 
	
		else 
			echo "Fastq for $samples_fastq_dir/$(basename $unmapped_bam .bam).fastq already exists - SKIP FASTQ CONVERSION (AND AUTOMATIC TRIMMING) "
		fi
		
				
		if [ ! -f $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed.fastq ]
		then
			echo "Start trimming $samples_fastq_dir/$(basename $unmapped_bam .bam).fastq"
			# trim fastq	
			
			# Trim first 6 base from each read, using Trimmomatic (HEADCROP:6)
			java -jar  $NGS_TRIMMOMATIC/trimmomatic-0.32-epignome.jar SE -phred33 -threads 30 -trimlog $samples_fastq_dir/r.log $samples_fastq_dir/$(basename $unmapped_bam .bam).fastq  $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed.fastq  HEADCROP:6 ILLUMINACLIP:$adapter_file:2:10:4:1:true:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16  
		
		else 
			echo "TRIMMED Fastq for $samples_fastq_dir/$(basename $unmapped_bam .bam).fastq already exists - SKIP FASTQ TRIMMING "
		fi
		
	else
	# PAIRED-END MODE
		# quick check to skip this step if fastq files exist already	
		if [ ! -f $samples_fastq_dir/$(basename $unmapped_bam .bam)_R1.fastq ]
		then
			echo "Start converting $unmapped_bam to $samples_fastq_dir/$(basename $unmapped_bam .bam)_R1.fastq"
		
			# Convert bam to fastq
# CB TEST: java -jar ~/tools/picard/SamToFastq.jar I=test1.bam F=test1.f.fastq F2=test1.r.fastq INCLUDE_NON_PF_READS=true
			java -jar $NGS_PICARD/SamToFastq.jar I=$unmapped_bam F=$samples_fastq_dir/$(basename $unmapped_bam .bam)_R1.fastq F2=$samples_fastq_dir/$(basename $unmapped_bam .bam)_R2.fastq INCLUDE_NON_PF_READS=true
	
			# Trim first 6 base from each read, using Trimmomatic (HEADCROP:6)
# CB TEST: java -jar ~/tools/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 -trimlog <samplename>_trimming.log -basein <samplename>.f.fastq -baseout <samplename>_trimmed.fastq HEADCROP:6 ILLUMINACLIP:epignome_adapters_2.fa:2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16

			java -jar  $NGS_TRIMMOMATIC/trimmomatic-0.32-epignome.jar PE -phred33 -threads 30 -trimlog $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimming.log $samples_fastq_dir/$(basename $unmapped_bam .bam)_R1.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_R2.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed_R1.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_unpaired_R1.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed_R2.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_unpaired_R2.fastq HEADCROP:6 ILLUMINACLIP:$adapter_file:2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16 
			#java -jar  $NGS_TRIMMOMATIC/trimmomatic-0.32.jar PE -phred33 -threads 30 -trimlog $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimming.log $samples_fastq_dir/$(basename $unmapped_bam .bam)_R1.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_R2.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed_R1.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_unpaired_R1.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed_R2.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_unpaired_R2.fastq HEADCROP:6 ILLUMINACLIP:$adapter_file:2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16 
		
		else 
			echo "Fastq for $samples_fastq_dir/$(basename $unmapped_bam .bam)_R1.fastq already exists  - SKIP FASTQ CONVERSION (AND AUTOMATIC TRIMMING)"
		fi
		
				# If fastq exists but not trimmed	
		if [ ! -f $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed_R1.fastq ]
		then
			echo "Start trimming $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed_R1.fastq"
		
			java -jar  $NGS_TRIMMOMATIC/trimmomatic-0.32-epignome.jar PE -phred33 -threads 30 -trimlog $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimming.log $samples_fastq_dir/$(basename $unmapped_bam .bam)_R1.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_R2.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed_R1.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_unpaired_R1.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed_R2.fastq $samples_fastq_dir/$(basename $unmapped_bam .bam)_unpaired_R2.fastq HEADCROP:6 ILLUMINACLIP:$adapter_file:2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16 
			
		else 
			echo "Fastq for $samples_fastq_dir/$(basename $unmapped_bam .bam)_R1.fastq already exists  - SKIP FASTQ TRIMMING"
		fi
		
	fi	


# -----------------------------------------------
#	BISMARK ALIGNMENT
# -----------------------------------------------

if [ 'ls -a $samples_fastq_dir' ]
then
	if [ $mode == "1" ]
	then
 fastq_file=$samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed.fastq 
	# Run Bismark on trimmed fastq files: alignment vs bisulfite genome and methylation calling
	
			bismark  $bismark_indexed_genome $fastq_file  --bam --unmapped --path_to_bowtie /cm/shared/apps/bowtie/2.2.3/bin --bowtie2 --temp_dir bismark_temp --output_dir $samples_aligned_dir -p 4 --prefix $genome
		

		
	else
fastq_file=$samples_fastq_dir/$(basename $unmapped_bam .bam)_trimmed_R1.fastq
	# PAIRED-END
			# Run Bismark on trimmed fastq files: alignment vs bisulfite genome and methylation calling
		
		# add --upto 1000000 for test run
			bismark  $bismark_indexed_genome --1 $fastq_file  --2 $samples_fastq_dir/$(basename $fastq_file _R1.fastq)_R2.fastq --bam --unmapped --path_to_bowtie /cm/shared/apps/bowtie/2.2.3/bin --bowtie2 --temp_dir bismark_temp --output_dir $samples_aligned_dir --minins 0 --maxins 5000 -p 4 --prefix $genome

# CB TEST: bismark /dseq/lab/bsf/scratch/anuzzo/meth-seq/bisulfite-genome/bowtie2-indexed/mm10/ --1 <samplename>_trimmed_1P.fastq --2 <samplename>_trimmed_2P.fastq <samplename>_trimmed_1U.fastq --bowtie2 --path_to_bowtie /slow/home/cbock/tools/bowtie2-2.1.0 --temp_dir Bismark_temp --minins 0 --maxins 5000 -p 8 --prefix <samplename>_bismark_mm10

			
	
	fi
fi

bismark_aligned_file=$genome.$(basename $unmapped_bam .bam)_trimmed.fastq_bismark_bt2.bam

if [ $mode == "2" ]
then
	bismark_aligned_file=$genome.$(basename $unmapped_bam .bam)_trimmed_R1.fastq_bismark_bt2_pe.bam
	
fi

#rename unmapped aligned $bismark_aligned_file

echo "bismark aligned filename: "
echo $bismark_aligned_file

#aligned_file=$(basename $unmapped_bam)
# -----------------------------------------------
#	REMOVE PCR DUPLICATE
# -----------------------------------------------

#for bam_file in `ls $samples_aligned_dir/*.bam`; do
if [ $mode == "1" ] 
	then
		deduplicate_bismark --single $samples_aligned_dir/$bismark_aligned_file --bam
	else
	# CB optim: deduplicate_bismark --paired test1_bismark_mm10.bam --bam --samtools_path ~/tools/samtools-0.1.19/
		deduplicate_bismark --paired $samples_aligned_dir/$bismark_aligned_file --bam 
	fi
	#samtools rmdup -S  $bam_file $samples_aligned_dir/nodup_$(basename $bam_file)
#done

echo 'PCR removal done.'


# -----------------------------------------------
#	METHYLATION CALL
# -----------------------------------------------

dedup_file=$(basename $bismark_aligned_file .bam).deduplicated.bam

#for bam_file in `ls $samples_aligned_dir/*deduplicated*.bam`; do
echo "dedup filename: "
echo $dedup_file
	if [ $mode == "1" ] 
	then
		bismark_methylation_extractor --single-end  --report --bedGraph --merge_non_CpG --cytosine_report --genome_folder $bismark_indexed_genome --gzip --output $meth_dir $samples_aligned_dir/$dedup_file > $meth_dir/meth-output.txt
	else
		bismark_methylation_extractor --paired-end --no_overlap --report --bedGraph --merge_non_CpG --cytosine_report --genome_folder $bismark_indexed_genome --gzip --output $meth_dir $samples_aligned_dir/$dedup_file > $meth_dir/meth-output.txt
	
	# convert bam file into sam file and sort again to compensate for a sorting issue of "deduplicate_bismark"
	fi
#dedup_file=$samples_fastq_dir/$(basename $bismark_aligned_file .bam).deduplicated.bam
samtools view $samples_aligned_dir/$dedup_file | LANG=C sort --temporary-directory=$project_folder/$samples_group/$pipeline_version/$genome/log/tmp  > $samples_aligned_dir/$(basename $dedup_file .bam).sam
	
#done



echo 'Methylation call done'

# -----------------------------------------------
#	ALIGNED READ FILTERING
# -----------------------------------------------

sam_file=$(basename $dedup_file .bam).sam

#for sam_file in `ls $samples_aligned_dir/*.sam`; do
echo "sam filename: "
echo $sam_file
	if [ $mode == "1" ] 
	then
		python $src_folder/bisulfiteReadFiltering.py --infile=$samples_aligned_dir/$sam_file --outfile=$samples_aligned_dir/$(basename $sam_file .sam).filtered.sam --genome=$genome --genomeDir=$ref_genome --minNonCpgSites=3 --minConversionRate=0.9
		echo 'filtering OK'
	else
			echo 'START: filtering'
		
		python $src_folder/bisulfiteReadFiltering.py --infile=$samples_aligned_dir/$sam_file --outfile=$samples_aligned_dir/$(basename $sam_file .sam).filtered.sam --genome=$genome --genomeDir=$ref_genome --pairedEnd --minNonCpgSites=3 --minConversionRate=0.9
	echo 'filtering OK'
	fi
#done

echo 'Filtering done' 
echo 'Sorting and indexing BAM done'

# -----------------------------------------------
#	FINAL SORTING AND INDEXING
# -----------------------------------------------
filtered_file=$(basename $sam_file .sam).filtered.sam

echo "filtered filename: "
echo $filtered_file

#for sam_file in `ls $samples_aligned_dir/*filtered.sam`; do
# create sorted and indexed BAM files for visualization and further analysis

	java -jar $NGS_PICARD/ReplaceSamHeader.jar I=$samples_aligned_dir/$filtered_file HEADER=$samples_aligned_dir/$(basename $filtered_file .filtered.sam).bam O=$samples_aligned_dir/$(basename $filtered_file .sam).bam
	java -jar $NGS_PICARD/SortSam.jar I=$samples_aligned_dir/$(basename $filtered_file .sam).bam O=$samples_aligned_dir/$(basename $filtered_file .sam)_final.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
	java -jar $NGS_PICARD/BuildBamIndex.jar I=$samples_aligned_dir/$(basename $filtered_file .sam)_final.bam 
	#samtools sort  $bam_file $samples_aligned_dir/$(basename $bam_file .bam)_sorted
	#samtools index $samples_aligned_dir/$(basename $bam_file .bam)_sorted.bam
#done




#-------------------------------
#  imported from wgbs-4.0.sh 
# -------------------------------

# output dir for bismark metilation extractor 
filtered_meth_dir=$samples_aligned_dir/extractor-post-filter
mkdir -p $filtered_meth_dir
chmod 775 $filtered_meth_dir

# Assumes bimsark alignement and duplicate removal have been done



filtered_bam=$samples_aligned_dir/$(basename $filtered_file .sam).bam
echo "filtered bam filename: "
echo $filtered_bam
	if [ $mode == "1" ] 
	then
		bismark_methylation_extractor --single-end  --report --bedGraph --merge_non_CpG --cytosine_report --genome_folder $bismark_indexed_genome --gzip --output $filtered_meth_dir $filtered_bam > $filtered_meth_dir/meth-output.txt
	else
		bismark_methylation_extractor --paired-end --no_overlap --report --bedGraph --merge_non_CpG --cytosine_report --genome_folder $bismark_indexed_genome  --gzip --output $filtered_meth_dir $filtered_bam > $filtered_meth_dir/meth-output.txt
	
	# convert bam file into sam file and sort again to compensate for a sorting issue of "deduplicate_bismark"
	fi
# dedup_file=$samples_fastq_dir/$(basename $bismark_aligned_file .bam).deduplicated.bam
# samtools view $samples_aligned_dir/$filtered_file | LANG=C sort > $samples_aligned_dir/$(basename $filtered_file .bam).sam

	# samtools view $filtered_file | LANG=C sort --temporary-directory=$project_folder/$samples_group/$pipeline_version/$genome/log/tmp > $samples_aligned_dir/$(basename $filtered_file .bam).sam
			
#done



echo 'POST-filter Methylation call done'



