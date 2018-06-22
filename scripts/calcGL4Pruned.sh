#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J calcGL
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcGL_%j_%A_%a.out
#SBATCH -t 24:00:00
#SBATCH --mem=30G
#SBATCH -n 4
#SBATCH --export=OUT_DIR=geno_lik/pass1/pruned_chunks/,IN_DIR=var_sites/pass1/pruned_positions/positions_,BAM_LIST=pass1_bam.all.list

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load angsd
module load angsd

# make directory to store output (if doesn't yet exist)
mkdir -p "$OUT_DIR"

# set variables
# pad the task id with leading zeros
printf -v TASK_ID "%03g" $SLURM_ARRAY_TASK_ID
# necessary becuase slurm array task ID doesn't hold leading zeros
POS_FILE=$IN_DIR"chunk"$TASK_ID
echo "position file: $POS_FILE"
echo "out directory: $OUT_DIR"

# apply filtering with SAMtools & PICARD
echo "finding genotype likelihood using ANGSD on BAMS for hilo pruned SNPs chunk $SLURM_ARRAY_TASK_ID"
# steps:
# (0) Start with filtered BAM files, reference genome and list of positions in each chunk
# (1) For each chunk, index positions with ANGSD
angsd sites index $POS_FILE
# (2) make a list of regions with positions in that chunk
# cut -f1 $POS_FILE | uniq > $POS_FILE.chr 
# above line wasn't specific enough b/c it reads whole chromosomes
# below specifies smaller regions, e.g. 1:2000-4500
rm -f $POS_FILE.chr
# first remove old file if it already exists
for i in $(cut -f1 $POS_FILE | uniq); \
do echo $i:$(awk -v i="$i" '$1 == i {print $2}' \
$POS_FILE | head -n 1)-$(awk -v i="$i" '$1 == i {print $2}' \
$POS_FILE | tail -n 1) >> $POS_FILE.chr; done
# (3) calculate genotype likelihoods using samtools algorithm and quality filters
echo "out: $OUT_DIR/chunk_$TASK_ID"

angsd -out "$OUT_DIR/chunk_$TASK_ID" \
-doMajorMinor 3 \
-sites $POS_FILE -rf $POS_FILE.chr \
-GL 1 -doGlf 2 \
-minMapQ 30 -minQ 20 \
-bam $BAM_LIST \
-remove_bads 1 \
-P 4

# settings:
# -sites and -rf specify which positions and chromosomes to calculate GL's for
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 3 takes major and minor allele from -sites file (for consistency)
# -bam list of bams to include (all newly sequenced allopatric mex. and sympatric mexicana & maize pops)
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -P 4 splits the analysis job over 4 nodes (but does not distribute I/O)
