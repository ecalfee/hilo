#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/genes_r_ancestry
#SBATCH -J defK
#SBATCH -o /home/ecalfee/hilo/slurm-log/defineKancestries_fromNGSadmix_output_%A_%a.out
#SBATCH -t 10:00
#SBATCH --mem=4G

# to run: sbatch --export=PREFIX=pass2_alloMAIZE_PalmarChico,WIND=1cM,K=3 --array=0-100 define_K_ancestries_from_NGSadmix.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load R

# run R script
echo "defining K ancestries"
for r in {1..5}
do Rscript ./define_K_ancestries.R ${PREFIX} ${r} ${WIND} ${BOOT} ${K}
done

echo 'all done!'
