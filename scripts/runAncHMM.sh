#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J ancHMM
#SBATCH -o /home/ecalfee/hilo/slurm-log/runAncHMMj_%A_%a.out
#SBATCH -t 192:00:00
#SBATCH --mem=24G
#SBATCH --array=366

# this script runs local ancestry inference using ancestry_hmm

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# directory with input files
DIR="var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/ancestry_hmm/input"

#LIST_OF_POPS=($(awk '{print $1}' landraces_fromLi/alloMaizeInclude.list)) # make array of individuals
#LIST_OF_ALPHA_MAIZE=
#LIST_OF_ALPHA_MEX=
#POP=${LIST_OF_INDS[$SLURM_ARRAY_TASK_ID]} # get individual i based on SLURM_ARRAY_TASK_ID
POP=pop${SLURM_ARRAY_TASK_ID} # e.g. pop366
ALPHA_MAIZE=0.86
ALPHA_MEX=0.14

# go to directory with input files
cd ${DIR}

#run ancestry_hmm
echo "running local ancestry inference pop"${POP}
ancestry_hmm -a 2 ${ALPHA_MAIZE} ${ALPHA_MEX} \
-p 0 100000 ${ALPHA_MAIZE} -p 1 -100 ${ALPHA_MEX} \
--ne 10000 --timin 0 --tmax 1000 \
-i ${POP}.anc_hmm.input \
-s ${POP}.anc_hmm.ids.ploidy

echo "all done!"
