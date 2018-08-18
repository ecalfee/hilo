#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J MAF2SitesFile
#SBATCH -o /home/ecalfee/hilo/slurm-log/MAF2SitesFile_%j_%A_%a.out
#SBATCH -t 30:00
#SBATCH --mem=200M
#SBATCH --array=0-425
#SBATCH --export=DIR=geno_lik/merged_pass2_all_alloMaize4Low_16/allVar

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset


echo "making sites file out of a mafs.gz file for region "$SLURM_ARRAY_TASK_ID
zcat ${DIR}/region_$SLURM_ARRAY_TASK_ID.mafs.gz | \
awk '$1 != "chromo" {print $1 "\t" $2 "\t" $3 "\t" $4}' \
${DIR}/region_$SLURM_ARRAY_TASK_ID.SNPs

echo "all done region "$SLURM_ARRAY_TASK_ID
