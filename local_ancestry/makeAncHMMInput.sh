#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/local_ancestry
#SBATCH -J inputHMM
#SBATCH -o /home/ecalfee/hilo/slurm-log/makeAncHMMInput_%A_%a.out
#SBATCH -t 30:00
#SBATCH --mem=4G
#SBATCH --array=18-19,21,23-31,34-35,360-363,365-374

# to run: sbatch --export=PREFIX=pass2_alloMAIZE makeAncHMMInput.sh

# array contains all the allopatric population #s
popN=$SLURM_ARRAY_TASK_ID

# paths to input and output directories
dir_allo_counts = "results/counts/$PREFIX"
dir_symp_counts = "results/countsMajMin/$PREFIX"
inds_file = "../global_ancestry/results/NGSAdmix/$PREFIX/globalAdmixtureByIncludedIndividual.txt" # has list of all included individuals, and which pop
dir_output = "results/ancestry_hmm/$PREFIX/input"


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
Rscript ../scripts/make_input_ancestry_hmm.R $popN $dir_allo_counts $dir_symp_counts $inds_file $dir_output

echo 'all done!'
