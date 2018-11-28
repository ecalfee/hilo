#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/scripts
#SBATCH -J windCM
#SBATCH -o /home/ecalfee/hilo/slurm-log/makeWindowsCM_%A_%a.out
#SBATCH -t 1:00:00
#SBATCH --mem=8G
#SBATCH --array=1-10

i=$SLURM_ARRAY_TASK_ID
CM_WINDOW=0.1
DIR_IN="../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM"
SNP_FILE_IN=$DIR_IN"/chr"$i".var.sites"
# produces output bed file in same directory, $DIR_IN"chr"$i""_"$CM_WINDOW"cM_windows.bed")

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load R

# run R script
echo "making cM windows around SNPs"
Rscript makeWindowsCM.R $CM_WINDOW $SNP_FILE_IN
echo 'all done!'
