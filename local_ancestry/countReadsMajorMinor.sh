#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J countMajMin
#SBATCH -o /home/ecalfee/hilo/slurm-log/countReadsMajorMinor_%A_%a.out
#SBATCH -t 2:00:00
#SBATCH --mem=4G
#SBATCH --array=1-200

# array contains all the chromosome #
ID=$SLURM_ARRAY_TASK_ID
# directory for input and output
dir="geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load R

# run R script for each chromosome
for CHR in {1..10}
    do echo "converting ACGT read counts to major/minor read counts for hilo"$ID" chr"$CHR
    Rscript ../scripts/count_reads_major_minor.R $CHR $ID $dir
  done
echo 'all done!'
