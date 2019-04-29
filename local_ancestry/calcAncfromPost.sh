#!/bin/bash
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/local_ancestry
#SBATCH -J calcAnc
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcAnc_%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --mem=8G
#SBATCH --array=18-19,21,23-31,34-35,360-363,365-374

# to run: sbatch --export=SUBDIR=output_noBoot,PREFIX=pass2_alloMAIZE calcAncfromPost.sh

# note: pop 366 is a good one to start with and has index 19

# this script takes posterior calls for all sites and all 3 genotypes
# and summarise it by population mexicana ancestry frequencies
# calc_genomewide_pop_anc_freq.R takes in HILOXX.posterior files from ancestry_hmm
# and returns and anc/ folder with individual alphas, ancestry frequencies,
# and pop ancestry frequencies
# Zanc_statistic.R is still early draft to run calculations from these files


# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load R
module load R

# directory with input/output subdirectories
DIR_MAIN="results/ancestry_hmm/$PREFIX/"
DIR_IN=${DIR_MAIN}/${SUBDIR}
DIR_OUT=${DIR_IN}/${anc}
# make output directory
mkdir -p ${DIR_OUT}

# use slurm array task id as the population number
POP_N=$SLURM_ARRAY_TASK_ID
IDs_FILE=${DIR_MAIN}/input/pop${POP_N}.anc_hmm.ids

#process output of ancestry_hmm
echo "summarising ancestry for pop"${POP_N}
Rscript calc_genomewide_pop_anc_freq.R ${POP_N} ${DIR_IN} ${IDs_FILE}

echo "all done!"
