#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J pruneFixBP
#SBATCH -o /home/ecalfee/hilo/slurm-log/pruneFixBP_%j_%A_%a.out
#SBATCH -t 16:00:00
#SBATCH --mem=8G
#SBATCH --array=1-10
#SBATCH --export=MIN_IND=75

# the default for MIN_IND is set to 75 but can be overwritten when calling the script
# e.g. sbatch --export=MIN_IND=100 pruneFixedBP.sh

# general bash script settings to make sure if any errors in the pipeline fail
set –o pipefail
set –o errexit
set –o nounset

# load modules
# using python3 for pandas
module load python3

# make output directory
mkdir -p var_sites/pass1/sympatric/pruned/minN$MIN_IND/

# prune to fixed distance of 10kb
python3 ../scripts/pruneFixedBP.py 10000 $MIN_IND 0.1 var_sites/pass1/sympatric/pruned/minN$MIN_IND/thin10kb_chr$SLURM_ARRAY_TASK_ID.mafs \
$(for i in $(awk -v chr=$SLURM_ARRAY_TASK_ID '$1 == chr {print $4}' refMaize/divide_50Mb/ALL_regions.list); do echo var_sites/pass1/sympatric/region_$i.mafs.gz; done
)
# second line finds all regions associated with a specific chromosome (from file ALL_regions.list) and lists those mafs.gz files as input files to pruneFixedBP.py
echo "done pruning SNPs by fixed BP for NGSadmix in hilo sympatric data only, minN="$MIN_IND" individuals, chr"$SLURM_ARRAY_TASK_ID
