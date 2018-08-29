#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J vcf2counts
#SBATCH -o /home/ecalfee/hilo/slurm-log/alloVCF2AlleleCounts_%j_%A_%a.out
#SBATCH -t 1:00:00
#SBATCH --mem=2G
#SBATCH --array=1-10

# array is the chromosome #
i=$SLURM_ARRAY_TASK_ID
POP="maize.allo.4Low16"
dir="var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/"

# this script takes in a VCF file (for an allopatric population)
# and outputs a counts file from plink

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load plink

echo "counting ref and alt allele in vcf chr"${i}" pop "${POP}
# take name.vcf.gz file from angsd and
# (1) convert it to a simple vcf format using plink, name.vcf (drops genotype likelihoods etc.)
# this file can then be read in by GATK for major/minor allele to get read counts for all low-coverage individuals
# (2) summarize allele counts of major and minor allele for high-coverage allopatric maize in a plink output file, name.frq.count
plink --vcf ${dir}/${POP}_chr${i}.vcf.gz --recode vcf-iid --freq counts --keep-allele-order --out ${dir}/${POP}_chr${i}

echo "all done!"
