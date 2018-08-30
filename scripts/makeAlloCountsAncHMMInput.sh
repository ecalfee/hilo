#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J alloHMM
#SBATCH -o /home/ecalfee/hilo/slurm-log/makeAlloCountsAncHMMInput_%A_%a.out
#SBATCH -t 30:00
#SBATCH --mem=4G
#SBATCH --array=1-10

# array contains all the chromosome #
CHR=$SLURM_ARRAY_TASK_ID

# directory for input and output
dir="../data/var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load R


# run R script
echo "making ancestry hmm input allopatric counts for chr"$CHR
Rscript ../scripts/make_allo_counts_ancestry_hmm.R $CHR $dir

echo 'all done!'
