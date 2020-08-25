#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/local_ancestry
#SBATCH -J vcf2counts
#SBATCH -o /home/ecalfee/hilo/slurm-log/alloVCF2AlleleCounts_%A_%a.out
#SBATCH -t 1:00:00
#SBATCH --mem=2G
#SBATCH --array=1-10

# to run: sbatch --export=PREFIX=pass2_alloMAIZE alloVCF2AlleleCounts.sh

# array is the chromosome #
i=$SLURM_ARRAY_TASK_ID
POP="allo.maize"
dir="results/counts/$PREFIX"

# this script takes in a VCF file (for an allopatric population)
# and outputs a counts file from plink

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load plink

echo "counting major and minor allele in vcf chr"${i}" pop "${POP}
# take name.vcf.gz file from angsd and
# (1) convert it to a simple vcf format using plink, name.vcf (drops genotype likelihoods etc.)
# (2) summarize allele counts of major and minor allele for high-coverage allopatric maize in a plink output file, name.frq.count
plink --vcf ${dir}/${POP}_chr${i}.vcf.gz --recode vcf-iid --freq counts --keep-allele-order --out ${dir}/${POP}_chr${i}

echo "all done!"
