#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/within_ancestry_diversity
#SBATCH -J addSites2Freq
#SBATCH -o /home/ecalfee/hilo/slurm-log/combinePopAlleleFreqWSites_%A_%a.out
#SBATCH -t 5:00
#SBATCH --mem=2G
#SBATCH -n 1
#SBATCH --array=0-425

# to run, e.g.
# sbatch --export="POP=parv,PREFIX=pass2_alloMAIZE" combinePopAlleleFreqWithSites.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

REGION=$SLURM_ARRAY_TASK_ID

# load R
module load R

echo "combining mafs.gz file with sites for pop $POP and region $REGION using snps from $PREFIX:"
Rscript ./combine_pop_allele_freqs_with_sites.R $POP results/allele_freq/$PREFIX/$POP/region_$REGION ../variant_sites/results/$PREFIX/region_$REGION.var.sites

echo "all done!"
