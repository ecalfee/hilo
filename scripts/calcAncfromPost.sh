#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/scripts
#SBATCH -J calcAnc
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcAnc_%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --mem=8G
#SBATCH --array=0-27
#SBATCH --export="SUBDIR=output_noBoot,GLOBAL_ADMIXTURE=globalAdmixtureByPopN.txt,ALL"

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
DIR_MAIN="../data/geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm"
DIR_IN=${DIR_MAIN}/${SUBDIR}
GLOBAL_ADMIXTURE_FILE="${DIR_MAIN}/input/${GLOBAL_ADMIXTURE}"
DIR_OUT=${DIR_IN}/${anc}
# make output directory
mkdir -p ${DIR_OUT}

# pull columns from file into arrays
LIST_OF_POPS=($(cut -d$'\t' -f 1  < $GLOBAL_ADMIXTURE_FILE))
LIST_OF_ALPHA_MAIZE=($(cut -d$'\t' -f 2  < $GLOBAL_ADMIXTURE_FILE))
LIST_OF_ALPHA_MEX=($(cut -d$'\t' -f 3  < $GLOBAL_ADMIXTURE_FILE))

# use slurm array task id as the array index to pull particular values out
i=$SLURM_ARRAY_TASK_ID
POP_N=${LIST_OF_POPS[${i}]} # e.g. 366 for pop366
ALPHA_MAIZE=${LIST_OF_ALPHA_MAIZE[${i}]}
ALPHA_MEX=${LIST_OF_ALPHA_MEX[${i}]}

# check for no admixture
echo "running pop"${POP_N}" maize: " ${ALPHA_MAIZE}" mex: "${ALPHA_MEX}
if [ ${ALPHA_MAIZE} = 0.00 ] || [ ${ALPHA_MAIZE} = 1.00 ]
then
    echo "pop"${POP_N}" is 100% one ancestry -> skipping local ancestry summary"
    exit 0 # exit without failure (no inference needed)
fi


#process output of ancestry_hmm
echo "summarising ancestry for pop"${POP_N}
Rscript calc_genomewide_pop_anc_freq.R ${POP_N} ${DIR_IN}
echo "all done!"
