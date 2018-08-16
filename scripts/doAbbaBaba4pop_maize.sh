#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J ABBA_BABA
#SBATCH -o /home/ecalfee/hilo/slurm-log/ABBA_BABA_%j_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=8G

# this script uses allopatric maize as the outgroup .. this will be underpowered to detect
# introgression between sympatric maize and mexicana if there is also gene flow to allopatric mexicana
OUT_FILE=AbbaBaba/4pop_maize_bam.regions.abbababa
BAM_FILE=AbbaBaba/4pop_maize_bam.list
SIZE_FILE=AbbaBaba/4pop_maize_bam.sizeFile
REGIONS=refMaize/random_regions/N1000.L100.regions
# to run: sbatch --export=POP1=18,POP2=19,DIR=SAF/pass1/N1000.L100.regions calc2DSFSAngsd.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

echo "now calculating ABBA-BABA"

angsd -doAbbababa2 1 \
-useLast 1 \
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
# -useLast 1 means make last bam file the 'outgroup'
# alternative without an outgroup is to do -Aanc which makes all ancestral alleles an A arbitrarily (and is otherwise an f3 test ?)
# -rf for regions file
# set basic mapQ and base Quality filters
# doing allele counts is required for abba-baba -doCounts 1
# block size I just made the same as my blocks in the genome..should be fine (?) .. they shouldn't have LD between them
