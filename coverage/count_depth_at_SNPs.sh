#!/bin/bash
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/coverage
#SBATCH -J countDepth
#SBATCH -o /home/ecalfee/hilo/slurm-log/countDepthAtSNPs_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=8G

# START WITH ZERO FOR ARRAY INDEXES (!)

# this script takes in a bam file and a bed file of 10kb windows across the genome
# and calculates # of reads from bam that have Q > 30 and overlap each windown in the bed file
# to run: sbatch --array=0 --export=BED_POS=../local_ancestry/results/thinnedSNPs/pass2_alloMAIZE/whole_genome.bed, DIR_OUT=results/thinnedSNPs/pass2_alloMAIZE,PREFIX=pass2_alloMAIZE count_depth_at_SNPs.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio3 # for bedtools and samtools

i=$SLURM_ARRAY_TASK_ID
bams=($(cat ${DIR_IN}/${PREFIX_LIST}_bams.list))
BAM=${bams["$i"]}
ids=($(cat ${DIR_IN}/${PREFIX_LIST}_IDs.list))
ID=${ids["$i"]}

# make output directory
mkdir -p ${DIR_OUT}

# calculate coverage using bedtools
echo "calculating coverage bam $BAM"
samtools view -b -q 30 "$BAM" | bedtools coverage -sorted -a "$BED_POS" -b stdin > "$DIR_OUT/$ID.bed"

echo "all done!"

# options
# samtools -q 30 limits to reads with mapping quality > 30
# -b specifies bam output
# bedtools coverage calculates number of reads in -b that overlap any part of features in -a
# -sorted assumes b is sorted by position to reduce memory cost (which these bams are)
