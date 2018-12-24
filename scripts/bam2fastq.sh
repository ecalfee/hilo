#!/bin/bash

# this script is a helper script called by slurm scripts
# with two input arguments: bam prefix and fastq prefix
# it takes bam files with reads aligned to one reference and converts them to fastq
# so these .fq files can then be remapped to the APGv4 reference
# to run: bam2fastq.sh input_dir/HILO23 output_dir/HILO23

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio # loads samtools and bedtools

echo "$0"
echo "$# parameters"; echo "$@"

BAM_IN_PREFIX=$1
FASTQ_OUT_PREFIX=$2

echo "sorting BAM by name"
# use samtools to sort input bam file by name (so paired reads are together)
samtools sort -n ${BAM_IN_PREFIX}".bam" ${BAM_IN_PREFIX}".sortName.bam"

echo "BAM -> fastq"
# use bedtools to convert bam to fastq file only if fastq does not already exist
bedtools bamtofastq -i ${BAM_IN_PREFIX}".sortName.bam" \
-fq ${FASTQ_OUT_PREFIX}"_1.fq" -fq2 ${FASTQ_OUT_PREFIX}"_2.fq"

echo "all done!"
# may want to delete intermediate file .sortName.bam later for space
