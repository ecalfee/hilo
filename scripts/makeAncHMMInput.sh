#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J inputHMM
#SBATCH -o /home/ecalfee/hilo/slurm-log/makeAncHMMInput_%A_%a.out
#SBATCH -t 30:00
#SBATCH --mem=4G
#SBATCH --array=18-19,21,23-31,34-35,360-363,365-374

# array contains all the allopatric population #s
popN=$SLURM_ARRAY_TASK_ID

input_dir="geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM"
output_dir="geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load R

# make directory to store output (if doesn't yet exist)
mkdir -p $output_dir

# run R script
echo "making ancestry hmm input files for population "$popN
Rscript ../scripts/make_input_ancestry_hmm.R $popN $input_dir $output_dir

echo 'all done!'
