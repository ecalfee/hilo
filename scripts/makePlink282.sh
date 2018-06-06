#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J plink282
#SBATCH -o /home/ecalfee/hilo/slurm-log/plink282_%j_%A_%a.out
#SBATCH -t 1:00:00
#SBATCH --mem=20G
#SBATCH --array=1-10

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load plink
module load plink

# make directory to store output (if doesn't yet exist)
mkdir -p var_sites/maize282/

# set variables
# pad the task id with leading zeros
printf -v TASK_ID "%03g" $SLURM_ARRAY_TASK_ID

# convert chromosomes 1-10 from vcf to plink bed/bim/fam for 282 diversity panel
# then filter for low LD
plink --vcf /group/jrigrp/Share/genotypes/282_7X/c"$SLURM_ARRAY_TASK_ID"_282_corrected_onHmp321.vcf.gz \
--out var_sites/maize282/chr$SLURM_ARRAY_TASK_ID \
--keep-allele-order --make-bed --biallelic-only list
# --keep-allele-order preserves the original ref/alt allele in the vcf 
# --biallelic-only list keeps only biallelic alleles and makes a list of skipped variant IDs to plink.skip.3allele
# example vcf file:
#/group/jrigrp/Share/genotypes/282_7X/c10_282_corrected_onHmp321.vcf.gz


