#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J ancNoBoot
#SBATCH -o /home/ecalfee/hilo/slurm-log/runAncHMM_noBoot_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH --mem=8G
#SBATCH --array=0-27
#SBATCH --export="Ne=10000,SUBDIR_OUT=output_noBoot,ALL"

# note: pop 366 is a good one to start with and has index 19
# try loading bio module
module load bio
module load Ancestry_HMM # loads copy from farm

# this script runs local ancestry inference using ancestry_hmm

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# directory with input/output subdirectories
DIR="geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm"
cd ${DIR} # move to main directory
GLOBAL_ADMIXTURE_FILE="input/globalAdmixtureByPopN.txt"

# pull columns from file into arrays
LIST_OF_POPS=($(cut -d$'\t' -f 1  < $GLOBAL_ADMIXTURE_FILE))
LIST_OF_ALPHA_MAIZE=($(cut -d$'\t' -f 2  < $GLOBAL_ADMIXTURE_FILE))
LIST_OF_ALPHA_MEX=($(cut -d$'\t' -f 3  < $GLOBAL_ADMIXTURE_FILE))

# use slurm array task id as the array index to pull particular values out
i=$SLURM_ARRAY_TASK_ID
POP=pop${LIST_OF_POPS[${i}]} # e.g. pop366
ALPHA_MAIZE=${LIST_OF_ALPHA_MAIZE[${i}]}
ALPHA_MEX=${LIST_OF_ALPHA_MEX[${i}]}

# check for no admixture
echo "running "${POP}" maize: " ${ALPHA_MAIZE}" mex: "${ALPHA_MEX}
if [ ${ALPHA_MAIZE} = 0.00 ] || [ ${ALPHA_MAIZE} = 1.00 ]
then
    echo ${POP}" is 100% one ancestry -> skipping local ancestry inference"
    exit 0 # exit without failure (no inference needed)
fi

# make and go to directory where ancestry_hmm should output files
mkdir -p ${SUBDIR_OUT}
cd ${SUBDIR_OUT} # change directory to output directory

#run ancestry_hmm
echo "running local ancestry inference "${POP}
ancestry_hmm -a 2 ${ALPHA_MAIZE} ${ALPHA_MEX} \
-p 0 100000 ${ALPHA_MAIZE} -p 1 -100 ${ALPHA_MEX} \
--ne ${Ne} --tmin 0 --tmax 10000 \
-i ../input/${POP}.anc_hmm.input \
-s ../input/${POP}.anc_hmm.ids.ploidy

echo "all done!"
