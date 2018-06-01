#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J pruneFixed
#SBATCH -o /home/ecalfee/hilo/slurm-log/pruneFixed_%j_%A_%a.out
#SBATCH -t 00:30:00
#SBATCH --mem=2G

# general bash script settings to make sure if any errors in the pipeline fail
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load python

# prune to fixed distance 10kb
python ../scripts/pruneFixed.py True 0 1 10000 var_sites/pass1/chr$SLURM_ARRAY_TASK_ID.mafs.gz var_sites/pass1/chr$SLURM_ARRAY_TASK_ID.pruned.mafs.gz

echo "done pruning SNPs for NGSadmix in hilo CHR$SLURM_ARRAY_TASK_ID"
