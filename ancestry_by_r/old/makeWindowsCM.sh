#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/genes_r_ancestry
#SBATCH -J windCM
#SBATCH -o /home/ecalfee/hilo/slurm-log/makeWindowsCM_%A_%a.out
#SBATCH -t 15:00
#SBATCH --mem=8G
#SBATCH --array=1-10

# to run: genes_r_ancestry$ sbatch --export=PREFIX=pass2_alloMAIZE,CM_WINDOW=0.1 makeWindowsCM.sh

i=$SLURM_ARRAY_TASK_ID
#CM_WINDOW=0.1
#DIR_IN="../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM"
DIR_IN="../local_ancestry/results/thinnedSNPs/"$PREFIX
SNP_FILE_IN=$DIR_IN"/chr"$i".var.sites"
DIR_OUT=$DIR_IN"/windows"$CM_WINDOW"cM"
BED_FILE_OUT=$DIR_OUT"/chr"$i".bed"
# produces output bed file in same directory, $DIR_IN"chr"$i""_"$CM_WINDOW"cM_windows.bed")

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load R

# make output DIRECTORY
mkdir -p $DIR_OUT

# run R script
echo "making cM windows around SNPs"
Rscript makeWindowsCM.R $CM_WINDOW $SNP_FILE_IN $BED_FILE_OUT
echo 'all done!'
