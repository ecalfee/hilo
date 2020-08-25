#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/local_ancestry
#SBATCH -J alloHMM
#SBATCH -o /home/ecalfee/hilo/slurm-log/makeAlloCountsAncHMMInput_%A_%a.out
#SBATCH -t 30:00
#SBATCH --mem=4G
#SBATCH --array=1-10

# to run: local_ancestry$ sbatch --export=PREFIX=pass2_alloMAIZE makeAlloCountsAncHMMInput.sh

# array contains all the chromosome #
CHR=$SLURM_ARRAY_TASK_ID
MEX_IDs="../samples/pass2_pops/allo.mexicana_IDs.list"
PREFIX="pass2_alloMAIZE"
DIR_MEX_COUNTS="results/countsMajMin/$PREFIX"
MAIZE_FILE="results/counts/$PREFIX/allo.maize_chr$CHR"
DIR_SITES="results/thinnedSNPs/$PREFIX"
DIR_OUT="results/counts/$PREFIX"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load R

# run R script
echo "making ancestry hmm input allopatric counts for chr"$CHR
Rscript ./make_allo_counts_ancestry_hmm.R $CHR $MEX_IDs $DIR_MEX_COUNTS $MAIZE_FILE $DIR_SITES $DIR_OUT

echo 'all done!'
