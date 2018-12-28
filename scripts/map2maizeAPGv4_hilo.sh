#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J map2v4
#SBATCH -o /home/ecalfee/hilo/slurm-log/map2maizeAPGv4HILO_%A_%a.out
#SBATCH -t 3-00:00:00
#SBATCH --mem=32G
#SBATCH -n 16
#SBATCH --array=1-200

# this script takes in fastq files with raw paired reads
# and aligns reads to APGv4 reference with all chromosome and contigs
# inferring the following variables based on HILO ID 1-200:
# sample ID, fastq file 1, fastq file 2, output bam directory
# to run: sbatch map2maizeAPGv4map2maizeAPGv4_hilo.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio # loads samtools and bwa

REF="refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
i=$SLURM_ARRAY_TASK_ID
ID=HILO${i}
FASTQ1=$(ls /group/jrigrp6/DanAlignments/${ID}/${ID}_*_1.fq.gz)
FASTQ2=$(ls /group/jrigrp6/DanAlignments/${ID}/${ID}_*_2.fq.gz)
OUTPUT_DIR="hilo_bam_mapped2v4_allContigs"

# make output directory
mkdir -p ${DIR_OUT}

# note: reference genome needs to already be indexed, e.g. bwa index ref.fa
echo "running bwa mem for sample "${ID}
# realign fastq using bwa to maize reference AGPv4 official release
bwa mem -t 16 -v 3 ${REF} ${FASTQ1} ${FASTQ2}  | \
samtools view -bS -o ${OUTPUT_DIR}/hilo_${i}.bam -

echo "all done!"
# options
# bwa mem:
# -t k specifies k threads
# -v 3 is for extra verbose (comments in addition to error messages)
# 2 input fastq files (for paired reads)
# samtools:
# - means read from standard in
# -bS means SAM in and BAM out
