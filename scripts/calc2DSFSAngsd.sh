#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J 2DSFS
#SBATCH -o /home/ecalfee/hilo/slurm-log/2DSFS_%j_%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --mem=2G
#SBATCH -n 6

# to run: sbatch --export=POP1=18,POP2=19 calc2DSFSAngsd.sh
# note: need to set POP1 and POP2 variables

# pops >=1000 are the larger groupings of maize sympatric, mexicana sympatric and mexicana allopatric inds
# run separately 1000, 2000, 3000

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make folder to save results (if doesn't already exist):
mkdir -p 2D_SFS/pass1

# calculate pairwise 2D SFS
realSFS SAF/pass1/pop$POP1/chrALL.saf.idx SAF/pass1/pop$POP2/chrALL.saf.idx -P 6 > 2D_SFS/pass1/pop$POP1.pop$POP2.saf.idx

echo "done calculating pairwise 2D SFS for pop"$POP1" and pop"$POP2
