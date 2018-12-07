#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/scripts
#SBATCH -J f4
#SBATCH -o /home/ecalfee/hilo/slurm-log/f4fromPopFreq_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH --mem=8G
#SBATCH --export="POP1=maize.allo.4Low16,POP2=maize.symp,POP3=mexicana.symp,POP4=mexicana.allo,ALL"


# directory for input and output
DIR="../data/geno_lik/merged_pass1_all_alloMaize4Low_16/allVar_depthFilt/popFreqs"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load R

# run R script
echo "calculating f4 statistic"
Rscript f4_from_pop_freq.R ${POP1} ${POP2} ${POP3} ${POP4} ${DIR}

echo 'all done!'
