#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/within_ancestry_diversity
#SBATCH -J piBackground
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcPiAngsdBackground_%A_%a.out
#SBATCH -t 24:00:00
#SBATCH -n 2
#SBATCH --mem=16G
#SBATCH --array=0-2,4-34

# this script calculates within-ancestry pi for a set of regions from a file
# to run:
# within_ancestry_diversity$ sbatch --export=REGIONS=N1000L10kb,PREFIX=mex2 piAngsd_background.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# loads angsd
module load bio

i=$SLURM_ARRAY_TASK_ID # 0 indexed!

REGIONS_FILE="../data/refMaize/random_regions/$REGIONS.regions"

POPS=(symp.mexicana symp.maize allo.mexicana allo.maize pop18 pop19 pop20 pop21 pop22 pop23 pop24 pop25 pop26 pop27 pop28 pop29 pop30 pop31 pop33 pop34 pop35 pop360 pop361 pop362 pop363 pop365 pop366 pop367 pop368 pop369 pop370 pop371 pop372 pop373 pop374)

POP="${POPS["$i"]}"

DIR_OUT="results/$REGIONS/$PREFIX"
REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"

echo "Calculating pi for regions: $REGIONS_FILE"
echo "POP: $POP"

# make directory to store output (if doesn't yet exist)
mkdir -p "$DIR_OUT"

echo "finding site allele frequencies"
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) Use all sites to estimate site allele frequency
angsd -out "$DIR_OUT/$POP" \
-anc "$REF" \
-fold 0 \
-underFlowProtect 1 \
-rf "$REGIONS_FILE" \
-bam "results/input/$PREFIX/pops/$POP.list" \
-remove_bads 1 -minMapQ 30 -minQ 20 \
-GL 1 \
-doSaf 1 \
-P 2

echo "done with SAF! calculating SFS"
realSFS "$DIR_OUT/$POP.saf.idx" -P 2 > "$DIR_OUT/$POP.sfs"

echo "done with SFS! calculating within-pop diversity 'thetas'"
angsd -out "$DIR_OUT/$POP" \
-anc "$REF" \
-doThetas 1 \
-doSaf 1 \
-pest "$DIR_OUT/$POP.sfs" \
-underFlowProtect 1 \
-rf "$REGIONS_FILE" \
-bam "results/input/$PREFIX/pops/$POP.list" \
-remove_bads 1 -minMapQ 30 -minQ 20 \
-GL 1 \
-P 2

echo "summarizing thetas for region and windows"
thetaStat do_stat "$DIR_OUT"/"$POP".thetas.idx -outnames "$DIR_OUT"/"$POP".thetasAll
thetaStat do_stat "$DIR_OUT"/"$POP".thetas.idx -win 5000 -step 1000 -outnames "$DIR_OUT"/"$POP".thetasWindows


# options
# basic quality filtering for reads
# anc polarizes SFS by reference genome
# underFlowProtect is necessary for large #s of bams
