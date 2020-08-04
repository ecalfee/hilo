#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/variant_sites
#SBATCH -J MAF2SitesFile
#SBATCH -o /home/ecalfee/hilo/slurm-log/MAF2SitesFile_%A_%a.out
#SBATCH -t 15:00
#SBATCH --mem=2G
#SBATCH --array=0-425

#to run: sbatch --export="DIR=results/hilo_alloMAIZE_MAIZE4LOW,ALL" mafsToSitesFile.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio

echo "making sites file out of a mafs.gz file for region "$SLURM_ARRAY_TASK_ID
zcat ${DIR}/region_$SLURM_ARRAY_TASK_ID.mafs.gz | \
awk '$1 != "chromo" {print $1 "\t" $2 "\t" $3 "\t" $4}' \
> ${DIR}/region_$SLURM_ARRAY_TASK_ID.var.sites

# index the file too so that ANGSD can use it as a sites file
echo "indexing: "${DIR}"/region_"$SLURM_ARRAY_TASK_ID".var.sites"
sleep 2s # wait 2 seconds before indexing so that index doesn't have same timestamp as sites file
angsd sites index ${DIR}/region_$SLURM_ARRAY_TASK_ID.var.sites
echo "all done region "$SLURM_ARRAY_TASK_ID
