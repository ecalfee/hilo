#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J catSAF
#SBATCH -o /home/ecalfee/hilo/slurm-log/catSAF_%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --mem=2G
#SBATCH --array=18-31,33-35,360-363,365-374,1000,2000,3000

# pops >=1000 are the larger groupings of maize sympatric, mexicana sympatric and mexicana allopatric inds
# run separately 1000, 2000, 3000

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# go to folder with results by chromosome
cd SAF/pass1/pop$SLURM_ARRAY_TASK_ID
# and concatenate results together
realSFS cat -outnames chrALL chr1.saf.idx chr2.saf.idx chr3.saf.idx chr4.saf.idx chr5.saf.idx chr6.saf.idx chr7.saf.idx chr8.saf.idx chr9.saf.idx chr10.saf.idx

echo "done merging SAF for pop "$SLURM_ARRAY_TASK_ID
