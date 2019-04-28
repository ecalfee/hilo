#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/local_ancestry
#SBATCH -J countMajMin
#SBATCH -o /home/ecalfee/hilo/slurm-log/countReadsMajorMinor_%A_%a.out
#SBATCH -t 2:00:00
#SBATCH --mem=4G

# to run: local_ancestry$ sbatch --array=0-217 --export=PREFIX=pass2_alloMAIZE countReadsMajorMinor.sh

i=$SLURM_ARRAY_TASK_ID # which hilo sample from list file
ids=($(cat ../samples/"$PREFIX"_IDs.list))
ID=${ids["$i"]} # sample ID
DIR_IN="results/countsACGT/$PREFIX/$ID"
DIR_OUT="results/countsMajMin/$PREFIX/$ID"
DIR_SITES="results/thinnedSNPs/$PREFIX/"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make output DIRECTORY
mkdir -p "$DIR_OUT"

# load software needed
module load R

# run R script for each chromosome
for CHR in {1..10}
    do echo "converting ACGT read counts to major/minor read counts for hilo"$ID" chr"$CHR
    Rscript ./count_reads_major_minor.R $CHR $ID $DIR_SITES $DIR_IN $DIR_OUT
  done
echo 'all done!'
