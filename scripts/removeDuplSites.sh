#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J rmDuplSites
#SBATCH -o /home/ecalfee/hilo/slurm-log/removeDuplSites_%j_%A_%a.out
#SBATCH -t 2:00:00
#SBATCH --mem=2G
#SBATCH -n 1
#SBATCH --array=0-425%8
#SBATCH --export=DIR_SITES=geno_lik/merged_pass1_all_alloMaize4Low_16/allVar,DIR_OUT=var_sites/merged_pass1_all_alloMaize4Low_16

# %k ensures only k jobs max run at one time

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

SITES_FILE=${DIR_SITES}/region_$SLURM_ARRAY_TASK_ID.var.sites
# load angsd -- don't load -- updated version 9.20 is local
#module load angsd

# make directory to store output (if doesn't yet exist)
mkdir -p $DIR_OUT

echo "removing second duplicate copy of sites from region "$SLURM_ARRAY_TASK_ID

# split into two files _a and _b based on half the total length of sites in the original sites file
split -a 1 -l $(($( wc -l "$SITES_FILE" | awk '{print $1}' ) / 2)) \
"$SITES_FILE" "${DIR_OUT}/region_${SLURM_ARRAY_TASK_ID}_"

mv "${DIR_OUT}/region_${SLURM_ARRAY_TASK_ID}_a" "${DIR_OUT}/region_${SLURM_ARRAY_TASK_ID}.var.sites"

# index new positions
sleep 2s # wait 2 seconds before indexing so that index doesn't have same timestamp as sites file
angsd sites index "${DIR_OUT}/region_${SLURM_ARRAY_TASK_ID}.var.sites"
