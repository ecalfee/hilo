#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J outSAFAngsd
#SBATCH -o /home/ecalfee/hilo/slurm-log/outlierSAFAngsd_%A_%a.out
#SBATCH -t 3:00:00
#SBATCH -n 1
#SBATCH --mem=2G

# note: must set the directory to print output and find regions.txt file
# and must set array as well when calling the script:
# to run:
# (where 8 is the number of lines in hap_groups.list)
# sbatch --array=1-8 --export="DIR=outliers/chr4/inv4m,ALL" calcSAFAngsd_outliers.sh
i=$(($SLURM_ARRAY_TASK_ID-1)) # because 0 indexed
LIST_OF_POPS=($(cat $DIR/hap_groups.list))
POP=$LIST_OF_POPS[$i]
LIST_OF_BAMS=($(cat $DIR/hap_groups.files))
BAM=$LIST_OF_BAMS[$i]
REGION_FILE="$DIR/regions.txt"
DIR_OUT="$DIR/SAF"
REF=refMaize/AGPv4_no_contigs/AGPv4.fa

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# loads angsd
module load bio

# make directory to store output (if doesn't yet exist)
mkdir -p $DIR_OUT

echo "finding site allele frequencies for pop $POP, listed at $BAM for outlier $DIR"
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) Use all sites to estimate site allele frequency
angsd -out "$DIR_OUT/$POP" \
-anc $REF \
-fold 0 \
-underFlowProtect 1 \
-doThetas 1 \
-rf $REGION_FILE \
-bam $BAM \
-remove_bads 1 -minMapQ 30 -minQ 20\
-GL 1 -dosaf 1 \
-P 1
echo "done with SAF!"

# options
# all the usual quality filtering for reads
# doThetas calculates thetas
# anc polarizes SFS by reference genome
# underFlowProtect is necessary for large #s of bams
