#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J ABBA_Aanc
#SBATCH -o /home/ecalfee/hilo/slurm-log/ABBA_BABA_Aanc_%j_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=8G
#SBATCH --export=PREFIX=3pop_alloMaize_sympMaize_alloMex_bam

# this script uses an arbitrary A as the ancestral allele/outgroup
OUT_FILE=AbbaBaba/$PREFIX.regions.abbababa
BAM_FILE=AbbaBaba/$PREFIX.list
SIZE_FILE=AbbaBaba/$PREFIX.sizeFile
REGIONS=refMaize/random_regions/N1000.L100.regions

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

echo "now calculating ABBA-BABA"

angsd -doAbbababa2 1 \
-Aanc 1 \
-rf $REGIONS \
-bam $BAM_FILE \
-sizeFile $SIZE_FILE \
-out $OUT_FILE \
-minQ 20 \
-minMapQ 30 \
-doCounts 1
#-blockSize 100 \

echo "done calculating ABBA-BABA"
# v2 uses all individuals
# -Aanc makes all ancestral alleles an A arbitrarily (and is otherwise an f3 test ?)
# -rf for regions file
# block size I just made the same as my blocks in the genome..should be fine (?) .. they shouldn't have LD between them