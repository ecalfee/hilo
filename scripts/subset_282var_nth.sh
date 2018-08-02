#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/maize282/vcf_AGPv4
#SBATCH -J sub282var
#SBATCH -o /home/ecalfee/hilo/slurm-log/subset282varNth_%j_%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --mem=8G
#SBATCH --export=N=1000
#SBATCH --array=1-10

# thins the variants in the 282 vcf to a small set with low LD appropriate for PCA
# and removes non-biallelic SNPs
module load vcftools

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

vcftools --gzvcf chr$SLURM_ARRAY_TASK_ID.vcf.gz --remove-indels --thin 10000 --recode --out subset_chr$SLURM_ARRAY_TASK_ID
# filter's variants so that they only include SNPs and no two variants within 10kb of each other
# --recode means output a new vcf with prefix --out
# can't filter for multiallelic variants (to remove)

echo "all done"
