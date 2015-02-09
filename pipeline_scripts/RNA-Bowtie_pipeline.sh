#! /bin/bash
#Author: Andreas Schoenegger <aschoenegger@cemm.oeaw.ac.at>
#Author: Johanna Klughammer <jklughammer@cemm.oeaw.ac.at>
#Version: 2014-10-20

#TODO: 
# - Read set is currently hardcoded to "all". Code changes are necessary to enable pf-only/npf-only reads.
# - Picard location is set in function picard_sam_to_fastq
# - Trim_galore location is set in function trim_galore
# - BSMAP location is set in function bsmap
# - samtools has to be in $PATH otherwise skipped
# - fastqc has to be in $PATH otherwise skipped

#------------------------------
#tools
#------------------------------
picard_path=/cm/shared/apps/picard-tools/1.118/
trim_galore=/cm/shared/apps/trim_galore/0.3.3/trim_galore
bowtie1=/fhgfs/groups/lab_bock/jklughammer/Tools/bowtie-1.1.1/bowtie
bowtie2=/cm/shared/apps/bowtie/2.2.3/bin/bowtie2
trimmomatic_path=/fhgfs/groups/lab_bock/shared/trimmomatic_debug


module load samtools/0.1.19
module load cutadapt/1.4.2
module load trim_galore/0.3.3



# -----------------------------------------
# Preparations
# -----------------------------------------

if [ $# -eq 0 ]; then
	echo "Params! 1: analysis_dir, 2: unmapped_bam"
	exit 1
fi

analysis_dir=$1
unmapped_bam=$2
ref_genome_fasta=$3
ref_genome_name=$4
scripts=$5
setup=$6
adapter_file=$scripts/adapters/epignome_adapters_2_add.fa
ERCC_fasta=$7
readLen=$8


echo "Genome: "
echo $ref_genome_name 

# Directory exists check
if [ ! -d $analysis_dir ]; then
	echo "Analysis directory $analysis_dir does not exist!. Exiting!"
	exit 1
fi

# File exists check
if [ ! -f $unmapped_bam ]; then
	echo "BAM file $unmapped_bam does not exist!. Exiting!"
	exit 1
fi


# -----------------------------------------
# picard_sam_to_fastq
# -----------------------------------------
function picard_sam_to_fastq {
	input_bam=$1
	output_fastq=$2
	setup=$3
	
	#picard_path=/fhgfs/groups/lab_bsf/meth-seq/src/tools/picard-tools-1.100
	#picard_path=/fast/opt/software/picard-tools-1.107
	if [ $setup = "PE" ]; then
	java -jar $picard_path/SamToFastq.jar I=$input_bam F=$output_fastq F2=${output_fastq/.fastq/_2nd.fastq}  INCLUDE_NON_PF_READS=true
	fi
	if [ $setup = "SE" ]; then
	java -jar $picard_path/SamToFastq.jar I=$input_bam F=$output_fastq INCLUDE_NON_PF_READS=true
	fi
}


# -----------------------------------------
# trim_galore
# -----------------------------------------
function run_trim_galore {
	input_fastq=$1
	output_dir=$2

	#trim_galore=/fhgfs/groups/lab_bsf/meth-seq/src/tools/trim_galore_v0.3.3/trim_galore
	#trim_galore=/fast/opt/ngs_tools/trim_galore-0.2.7/trim_galore
	

	#trimGalore parameters
	#Adapter
	a="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	#quality trimming
	q=20
	#stringency: Overlap with adapter sequence required to trim a sequence
	s=1
	#Maximum allowed error rate
	e=0.1
	#Minimum Read length
	#by unchangeable default Trimmomatic discards reads of lenth 0 (produced by ILLUMINACLIP):
	#in order to be comparable to trimgalore MINLEN and --length were set to 16 in the respective tools
	l=16
	echo $q $s $e $l

	#   --- $trim_galore -q $q --phred33 -a $a --stringency $s -e $e --length $l --output_dir $output_dir $input_fastq 
	trim_galore -q $q --phred33 -a $a --stringency $s -e $e --length $l --output_dir $output_dir $input_fastq 
	

# --rrbs : removed after 05.03.2014 documentation (Optimization_conclusion.txt)
}

#-------------------------------------------
#trimmomatic
#-------------------------------------------
function run_trimmomatic {
	input_fastq1=$1
	input_fastq2=$2
	trimmed_fastq1=$3
	trimmedt_fastq2=$3
	output_dir=$(dirname $1)
	adapter_file=$5


	java -jar  $trimmomatic_path/trimmomatic-0.32-epignome.jar PE -phred33 -threads 6 -trimlog $output_dir/$(basename $input_fastq1 .fastq)_trimming.log $input_fastq1 $input_fastq2 $trimmed_fastq1 ${trimmed_fastq1/trimmed/unpaired} $trimmed_fastq2 ${trimmed_fastq2/trimmed/unpaired} HEADCROP:6 ILLUMINACLIP:$adapter_file:2:10:4:1:false:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16
	#java -jar  $trimmomatic_path/trimmomatic-0.32-epignome.jar PE -phred33 -threads 6 -trimlog $output_dir/$(basename $input_fastq1 .fastq)_trimming.log $input_fastq1 $input_fastq2 $trimmed_fastq1 ${trimmed_fastq1/trimmed/unpaired} $trimmed_fastq2 ${trimmed_fastq2/trimmed/unpaired} HEADCROP:6 ILLUMINACLIP:$adapter_file:2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:16

}

# -----------------------------------------
# contamination_check
# -----------------------------------------
function contamination_check {	
	#count FASTQ reads containing the adapter
	in_fastq=$1
	contamin_reads=`grep AGATCGGAAGAGC $in_fastq | wc -l`
	fastq_lines=`cat $in_fastq | wc -l`
	total_reads=`echo $(($fastq_lines/4))`
	#percent=`printf "%.2f" $(echo "$contamin_reads/$total_reads*100" | bc -l)`	# doesn't work if locale is DE.. (,/. issues)
	percent=`echo $contamin_reads $total_reads | awk '{ printf("%.2f", $1/$2*100) }'`
	
	echo "Contamination:$contamin_reads:$total_reads:$percent%" | tee -a $in_fastq.adapter.log
}


# -----------------------------------------
# Bowtie1
# -----------------------------------------
function run_bowtie {
	input_fastq=$1
	output_sam=$2
	ref_fasta=$3
	setup=$4

	echo "BOWTIE Input: $input_fastq"
	echo "BOWTIE outout: $output_bam"
	echo "BOWTIE ref genome: $ref_genome_fasta"
		
	#command="$bowtie2 -p 4 -L 15 -k 100 --local -x $ref_genome_fasta $input_fastq $output_sam"
	if [ $setup = "SE" ]; then
	fq=${input_fastq}_trimmed.fq
	command="$bowtie1 -q -p 6 -a -m 100 --sam $ref_fasta $fq $output_sam"
	fi
	
	if [ $setup = "PE" ]; then
	fq1=${input_fastq}_tc_trimmed.fq
	fq2=${input_fastq}_2nd_tc_trimmed.fq
	
	command="$bowtie1 -q -p 6 -a -m 100 --minins 0 --maxins 5000 --fr --sam $ref_fasta -1 $fq1 -2 $fq2 $output_sam"
	fi
	
	echo "COMMAND TO BE RUN: "
	echo $command
	
	`$command`

	#samtools flagstat $output_bam > $output_bam.flagstat
}





# -----------------------------------------
# Main
# -----------------------------------------

# --- 0a. Determine sample name ---
# sample=$(basename $unmapped_bam | cut -f2 -d# | cut -f1 -d.)
sample=$(basename $unmapped_bam | cut -f1 -d.)

echo $sample

# 0b. Create sample directory
fastq_dir=$analysis_dir/fastq/$sample
mkdir -p $fastq_dir

# --- 0c. QC ---
#>>>>>>>  samtools flagstat $unmapped_bam > $unmapped_bam.flagstat
# fastqc $unmapped_bam

# --- 1. Picard - SamToFastq ---
fastq=$fastq_dir/$sample.all.untrimmed.fastq
echo "Fastq file: $fastq"
if [ ! -f $fastq ]; then
	picard_sam_to_fastq $unmapped_bam $fastq $setup
fi

# --- 2. Trim Galore ---
trim_galore_outdir=$(dirname $fastq)

if [ $setup = "SE" ]; then
trimmed_fastq=${fastq/.fastq/_trimmed.fq} # Note: this path is set by trim_galore!
	if [ ! -f $trimmed_fastq ]; then
		contamination_check $fastq
		run_trim_galore $fastq $trim_galore_outdir
		contamination_check $trimmed_fastq
	fi
fi

if [ $setup = "PE" ]; then
	trimmed_fastq1=${fastq/.fastq/_tc_trimmed.fq}
	trimmed_fastq2=${fastq/.fastq/_2nd_tc_trimmed.fq}
	fastq2=${fastq/.fastq/_2nd.fastq}
	if [ ! -f $trimmed_fastq1 ]; then
		contamination_check $fastq
		contamination_check $fastq2
		run_trimmomatic $fastq $fastq2 $trimmed_fastq1 $trimmed_fastq2 $adapter_file
		contamination_check $trimmed_fastq1
		contamination_check $trimmed_fastq2
	fi
fi

## --- 3.1 Alignment with Tophat: ---

trimmed_fq_base=$(dirname $fastq)/$(basename $fastq .fastq)
logdir=$(dirname $(dirname $fastq))/log

sbatch --export=NONE --get-user-env=L --job-name=RNA-Th_$sample --ntasks=1 --cpus-per-task=2 --mem-per-cpu=8000 --partition=longq --time=3-00:00:00 -e $logdir/RNA_Th_${setup}_${sample}_%j.err -o $logdir/RNA_Th_${setup}_${sample}_%j.log $scripts/RNA-Tophat_pipeline.sh $trimmed_fq_base $ref_genome_name $setup $readLen


# --- 3.2 Alignment with Bowtie1: ---
align_dir="$analysis_dir/${ref_genome_name}_$setup/$sample"
mkdir -p $align_dir

aligned_sam=$align_dir/${sample}_bt1.sam
#>>>>> 

run_bowtie $trimmed_fq_base $aligned_sam $ref_genome_fasta $setup


#---4. make bam and flagstat
echo $(basename $aligned_sam .sam).bam
samtools view -bS $aligned_sam > ${aligned_sam/sam/"bam"}
samtools flagstat ${aligned_sam/sam/"bam"}  > ${aligned_sam/sam/"stats"}
samtools sort ${aligned_sam/sam/"bam"} ${aligned_sam/.sam/"_sorted"}
samtools index ${aligned_sam/.sam/"_sorted.bam"}
samtools depth ${aligned_sam/.sam/"_sorted.bam"} > ${aligned_sam/.sam/"_sorted.depth"}
#---4.1. filter reads to remove reads originating from genomic DNA

headerLines=`samtools view -SH $aligned_sam|wc -l`
if [ $setup = "SE" ]; then
python /fhgfs/groups/lab_bock/jklughammer/pipelines/bisulfiteReadFiltering_forRNA.py --infile=$aligned_sam --outfile=${aligned_sam/.sam/_filtered.sam} --skipped=${aligned_sam/.sam/_skipped.sam} --skipHeaderLines=$headerLines --minNonCpgSites=3 --minConversionRate=0.9 --maxConversionRate=0.1 --genome=$ref_genome_name --genomeDir=/fhgfs/groups/lab_bock/jklughammer/projects/comparative_epigenomics/genomes/ -r
fi

if [ $setup = "PE" ]; then
python /fhgfs/groups/lab_bock/jklughammer/pipelines/bisulfiteReadFiltering_forRNA.py --infile=$aligned_sam --outfile=${aligned_sam/.sam/_filtered.sam} --skipped=${aligned_sam/.sam/_skipped.sam} --skipHeaderLines=$headerLines --minNonCpgSites=3 --minConversionRate=0.9 --maxConversionRate=0.1 --pairedEnd --genome=$ref_genome_name --genomeDir=/fhgfs/groups/lab_bock/jklughammer/projects/comparative_epigenomics/genomes/ -r
fi

#change flags in skipped to unmapped
filtered_sam=${aligned_sam/.sam/_filtFlag.sam}
samtools view -hS ${aligned_sam/.sam/_filtered.sam} >$filtered_sam
awk -v skip=$headerLines -v OFS="\t" '{if (NR>skip){$2=4; print}}' ${aligned_sam/.sam/_skipped.sam} >> $filtered_sam



#---5. map unmapped reads to ERCCs
if [ $ERCC_fasta != "" ]; then
unmappable_bam=${aligned_sam/.sam/"_unmappable.bam"}
samtools view -b -f4 ${aligned_sam/sam/"bam"} > $unmappable_bam
if [ $setup = "SE" ]; then
fastq=${unmappable_bam/.bam/"_trimmed.fq"}
picard_sam_to_fastq $unmappable_bam $fastq $setup
trimmed_fq_base=$(dirname $fastq)/$(basename $fastq _trimmed.fq)
fi

if [ $setup = "PE" ]; then
fastq=${unmappable_bam/.bam/".fastq"}
picard_sam_to_fastq $unmappable_bam $fastq $setup
mv $fastq ${fastq/.fastq/_tc_trimmed.fq}
mv ${fastq/.fastq/_2nd.fastq} ${fastq/.fastq/_2nd_tc_trimmed.fq}
trimmed_fq_base=$(dirname $fastq)/$(basename $fastq .fastq)
fi


ERCC_sam=$align_dir/${sample}_bt1_ERCC.sam
run_bowtie $trimmed_fq_base $ERCC_sam $ERCC_fasta $setup
samtools view -bS $ERCC_sam|samtools sort - ${ERCC_sam/.sam/"_sorted"}
samtools index  ${ERCC_sam/.sam/"_sorted.bam"}


output_dir=$(dirname $aligned_sam)/bitSeq_bt1
mkdir -p $output_dir
Rscript $scripts/bitSeq_parallel.R $ERCC_sam $(dirname $(dirname $ERCC_fasta))/$(basename $ERCC_fasta).fa

fi


#---6. bitseq
output_dir=$(dirname $aligned_sam)/bitSeq_bt1
mkdir -p $output_dir

#Rscript $scripts/bitSeq_parallel.R $filtered_sam $(dirname $(dirname $ref_genome_fasta))/$ref_genome_name.fa
Rscript $scripts/bitSeq_parallel.R $aligned_sam $(dirname $(dirname $ref_genome_fasta))/$ref_genome_name.fa


