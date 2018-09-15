#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J 2DSFS
#SBATCH -o /home/ecalfee/hilo/slurm-log/2DSFS_%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --mem=1G
#SBATCH -n 6

# to run: sbatch --export=POP1=18,POP2=19,DIR=SAF/pass1/N1000.L100.regions calc2DSFSAngsd.sh
# note: need to set POP1 and POP2 and directory variables

# pops >=1000 are the larger groupings of maize sympatric, mexicana sympatric and mexicana allopatric inds
# run separately 1000, 2000, 3000

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# calculate pairwise 2D SFS
realSFS $DIR/pop$POP1.saf.idx $DIR/pop$POP2.saf.idx -P 6 > $DIR/pop$POP1.pop$POP2.sfs

echo "done calculating pairwise 2D SFS for pop"$POP1" and pop"$POP2
