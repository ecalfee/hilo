#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/scripts/
#SBATCH -J pruneFixcM
#SBATCH -o /home/ecalfee/hilo/slurm-log/pruneFixcM_%j_%A_%a.out
#SBATCH -t 16:00:00
#SBATCH --mem=8G
#SBATCH --array=75,100

# general bash script settings to make sure if any errors in the pipeline fail
set –o pipefail
set –o errexit
set –o nounset

# load modules
# using python3 because anaconda distribution has pandas
module load anaconda3

# prune to fixed distance .01cM (~10kb)
python pruneFixedcM.py 0.01 $SLURM_ARRAY_TASK_ID 0.1 ../data/var_sites/pass1/sympatric/pruned01cMMAF01minInd$SLURM_ARRAY_TASK_ID.mafs \
$(for i in {0..46}; do echo ../data/var_sites/pass1/sympatric/region_$i.mafs.gz; done)
echo "done pruning SNPs for NGSadmix in hilo sympatric data only"
