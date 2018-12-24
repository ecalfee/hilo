#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J map2v4
#SBATCH -o /home/ecalfee/hilo/slurm-log/map2maizeAPGv4BWA_%A_%a.out
#SBATCH -t 3-00:00:00
#SBATCH --mem=32G
#SBATCH -n 16

# this helper script takes in fastq files with raw paired reads
# and aligns reads to APGv4 reference with all chromosome and contigs
# designed to be called from a slurm script with two input arguments:
# sample ID, fastq file 1, fastq file 2, output bam directory
# to run: sbatch --export="ID=HILO23,FASTQ1=fastq_dir/HILO23_1.fq,FASTQ2=fastq_dir/HILO23_2.fq,DIR_OUT=../data/filtered_bams,ALL" map2maizeAPGv4.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio # loads samtools and bwa

REF="refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"

# make output directory
mkdir -p ${DIR_OUT}

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
