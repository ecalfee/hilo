#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/scripts
#SBATCH -J pipe
#SBATCH -o /home/ecalfee/hilo/slurm-log/alleleFreqPipeline_%j_%A_%a.out
#SBATCH -t 15:00
#SBATCH --mem=200M
#SBATCH --array=0-31
# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

POPS=(mexicana.allo.withXochi35 maize.allo.4Low16 pop18 pop19 pop20 pop21 pop22 pop23 pop24 pop25 pop26 pop27 pop28 pop29 pop30 pop31 pop34 pop35 pop360 pop361 pop362 pop363 pop365 pop366 pop367 pop368 pop369 pop370 pop371 pop372 pop373 pop374)

# to run
# sbatch calcAlleleFreqPop.sh pop22
# sbatch calcAlleleFreqPop.sh maize.symp
echo "now running for "${POPS[$SLURM_ARRAY_TASK_ID]}
sbatch --export=POP=${POPS[$SLURM_ARRAY_TASK_ID]},DIR_POPS=pass1_bam_pops,DIR_REGIONS=refMaize/divide_5Mb,DIR_SITES=geno_lik/merged_pass1_all_alloMaize4Low_16/allVar calcAlleleFreqPop.sh

echo "all done!"
