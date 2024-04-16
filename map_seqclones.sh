#!/bin/bash
#SBATCH --job-name=WGS_map_reads
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32GB
#SBATCH --array=0-95   ## we have 96 fastq files to process from two 96-plex runs 


##set variables
## this reference contains the GFP kanmx GAP1 CNV reporter 
REF='/scratch/work/cgsb/genomes/In_house/Fungi/Saccharomyces_cerevisiae/Gresham/GCF_000146045.2_R64_GAP1/GCF_000146045.2_R64_genomic_GAP1.fna'

##capture the output of a command line and store it in a variable
## run 1  - 97+97 fastqs
R1=($(ls /scratch/cgsb/gencore/out/Gresham/2022-08-15_HKCGLBGXM/merged/HKCGLBGXM_n01_*)) ### find a way for it to ignore the undetermined
R2=($(ls /scratch/cgsb/gencore/out/Gresham/2022-08-15_HKCGLBGXM/merged/HKCGLBGXM_n02_*))

##for the second run
##  96 + 96 fastqs. 
RUN2_R1=($(ls /scratch/cgsb/gencore/out/Gresham/2022-09-15_HGG2JDRX2/merged/HGG2JDRX2_n01_DGY*))
RUN2_R2=($(ls /scratch/cgsb/gencore/out/Gresham/2022-09-15_HGG2JDRX2/merged/HGG2JDRX2_n02_DGY*))

##this is going to assign the variables to file names
FORWARD1=${R1[$SLURM_ARRAY_TASK_ID]}
REVERSE1=${R2[$SLURM_ARRAY_TASK_ID]}

FORWARD2=${RUN2_R1[$SLURM_ARRAY_TASK_ID]}
REVERSE2=${RUN2_R2[$SLURM_ARRAY_TASK_ID]}

NAME=${FORWARD1:76:-9}

mkdir d_${NAME}
cd d_${NAME}

## trim adaptor sequences
module load cutadapt/3.1
## here are 4 adapter sequences from nextera and truseq protocols in the order: nextera/transposon seq read 2 - truseq read 1 - nextera read 1/transposon sequence - truseq adapter read 2
cutadapt -e 0.12 -m 30 -a CTGTCTCTTATACACATCT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATGTGTATAAGAGACAG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o trimmed_${NAME}_run1_R1.fastq -p trimmed_${NAME}_run1_R2.fastq ${FORWARD1} ${REVERSE1}
cutadapt -e 0.12 -m 30 -a CTGTCTCTTATACACATCT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATGTGTATAAGAGACAG -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o trimmed_${NAME}_run2_R1.fastq -p trimmed_${NAME}_run2_R2.fastq ${FORWARD2} ${REVERSE2}


##first, align reads with bwa mem
module load bwa/intel/0.7.17
bwa mem $REF trimmed_${NAME}_run1_R1.fastq  trimmed_${NAME}_run1_R2.fastq > ${NAME}_run1.bam

bwa mem $REF trimmed_${NAME}_run2_R1.fastq  trimmed_${NAME}_run2_R2.fastq > ${NAME}_run2.bam


##sort, index
module load samtools/intel/1.14
samtools sort ${NAME}_run1.bam > ${NAME}_run1_sort.bam
samtools index ${NAME}_run1_sort.bam 

samtools sort ${NAME}_run2.bam > ${NAME}_run2_sort.bam
samtools index ${NAME}_run2_sort.bam 

## merge bams
samtools merge -o ${NAME}_merged_sort.bam ${NAME}_run1_sort.bam ${NAME}_run2_sort.bam 

## obtaining alignment metrics using Picards tools
module purge
module load picard/2.23.8

java -jar /share/apps/picard/2.23.8/picard.jar MarkDuplicates I=${NAME}_merged_sort.bam O=${NAME}_marked_duplicates.bam M=${NAME}_marked_dup_metrics.txt

java -jar /share/apps/picard/2.23.8/picard.jar CollectAlignmentSummaryMetrics R=$REF I=${NAME}_marked_duplicates.bam O=${NAME}_alignment_metrics.txt 

## obtaining insert size metrics using Picards tools
java -jar /share/apps/picard/2.23.8/picard.jar CollectInsertSizeMetrics INPUT=${NAME}_marked_duplicates.bam OUTPUT=${NAME}_insert_metrics.txt HISTOGRAM_FILE=${NAME}_insert_size_histogram.pdf 

## obtaining read depth ie coverage using samtools
module load samtools/intel/1.14
samtools coverage ${NAME}_marked_duplicates.bam -o ${NAME}_coverage.txt
samtools depth -a ${NAME}_marked_duplicates.bam > ${NAME}_RD.txt


