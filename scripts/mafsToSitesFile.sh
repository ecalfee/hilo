#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J MAF2SitesFile
#SBATCH -o /home/ecalfee/hilo/slurm-log/MAF2SitesFile_%A_%a.out
#SBATCH -t 1:00:00
#SBATCH --mem=250M
#SBATCH --array=0-425
#SBATCH --export=DIR=geno_lik/merged_pass1_all_alloMaize4Low_16/allVar_depthFilt

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset


echo "making sites file out of a mafs.gz file for region "$SLURM_ARRAY_TASK_ID
zcat ${DIR}/region_$SLURM_ARRAY_TASK_ID.mafs.gz | \
awk '$1 != "chromo" {print $1 "\t" $2 "\t" $3 "\t" $4}' \
> ${DIR}/region_$SLURM_ARRAY_TASK_ID.var.sites

# index the file too so that ANGSD can use it as a sits file
echo "indexing"${DIR}"/region_"$SLURM_ARRAY_TASK_ID".var.sites"
sleep 2s # wait 2 seconds before indexing so that index doesn't have same timestamp as sites file
angsd sites index ${DIR}/region_$SLURM_ARRAY_TASK_ID.var.sites
echo "all done region "$SLURM_ARRAY_TASK_ID
