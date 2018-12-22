#!/bin/bash

# this helper script takes in fastq files with raw paired reads
# and aligns reads to APGv4 reference with all chromosome and contigs
# designed to be called from a slurm script with two input arguments:
# sample ID, fastq file 1, fastq file 2, output bam directory
# to run: map2maizeAPGv4.sh HILO23 fastq_dir/HILO23_1.fq fastq_dir/HILO23_2.fq ../data/filtered_bams

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio # loads samtools and bwa

ID=$1
FASTQ1=$2
FASTQ2=$3
OUTPUT_DIR=$4
REF="refMaize/AGPv4_official.fa"

# make output directory
mkdir -p $OUTPUT_DIR


# note: reference genome needs to already be indexed, e.g. bwa index ref.fa
echo "running bwa mem for sample "${ID}
# realign fastq using bwa to maize reference AGPv4 official release
bwa mem -t 16 -v 3 ${REF} ${FASTQ1} ${FASTQ2}  | \
samtools view -bS -o ${OUTPUT_DIR}/${ID}.bam -

echo "all done!"
# options
# bwa mem:
# -t k specifies k threads
# -v 3 is for extra verbose (comments in addition to error messages)
# 2 input fastq files (for paired reads)
# samtools:
# - means read from standard in
# -bS means SAM in and BAM out
