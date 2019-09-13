#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/filtered_bams
#SBATCH -J bam2fq
#SBATCH -o /home/ecalfee/hilo/slurm-log/bam2fastq_batch_%A_%a.out
#SBATCH -t 8-00:00:00
#SBATCH --mem=32G
#SBATCH -n 16

# this script takes in a bam file and creates 3 fastq files --
# PREFIX_1.fq.gz PREFIX_2.fq.gz PREFIX_3.fg.gz
# where _1 and _2 are the paired reads and _3 is the set of singletons
# before creating the fastq files, reads must be sorted by query name so pairs are in the same order

# to run: sbatch --array=0-13 --export=ID_file=../data/landraces_fromLi/alloMAIZE_IDs.list,DIR_IN=../data/landraces_fromLi/original,BAM=RIMMA0625_Andean.IndelRealigned.bam bam2fastq2v4.sh

i=$SLURM_ARRAY_TASK_ID
ids=($(cat ${DIR_IN}/${ID_file}))
ID=${ids["$i"]}

# so these .fq files can then be remapped to the APGv4 reference
mkdir -p "$DIR_IN"/fastq


# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio # loads samtools and bedtools

echo "sorting BAM by name then output fastq"
# use samtools to sort input bam file by name (so paired reads are together)
# then use samtools to output fastq files
samtools sort -n "$DIR_IN"/"$BAM" | \
samtools fastq -c 6 -1 "$DIR_IN"/fastq/"$ID"_1.fq.gz -2 "$DIR_IN"/fastq/"$ID"_2.fq.gz -s "$DIR_IN"/fastq/"$ID"_3.fq.gz -


echo "all done!"

# options

# samtools sort
# -n is for sorting by read name

# samtools fastq
# -s is for singletons, -1 for read 1 and -2 for read 2
# -c is for compression, level 6 is the default for gzip
