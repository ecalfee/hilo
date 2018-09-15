#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J GL282sites
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcGL4maize282sites_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH --mem=16G
#SBATCH -n 2
#SBATCH --array=1-10

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load angsd -- don't load angsd (local copy is most up-to-date v9.20
#module load angsd

# set variables
CHR="$SLURM_ARRAY_TASK_ID"
PREFIX="subset_chr$SLURM_ARRAY_TASK_ID"
SITES_DIR="maize282/vcf_AGPv4"
OUT_DIR="geno_lik/merged_pass1_all_alloMaize4Low_16"
BAM_LIST="merged_bam.pass1_all.alloMaize4Low_16.list"

echo "position file: $SITES_DIR/$PREFIX.sites"
echo "out directory: $OUT_DIR"

# make directory to store output (if doesn't yet exist)
mkdir -p "$OUT_DIR"

# calculate genotype likelihoods
# for minor and major allele specified in sites file (-doMajorMinor 3)
# -rf just speeds things up by using bam indexing
# first index sites file
angsd sites index "$SITES_DIR/$PREFIX.sites"
# then run angsd
angsd -out "$OUT_DIR/$PREFIX" \
-doMajorMinor 3 \
-minInd 0 \
-sites "$SITES_DIR/$PREFIX.sites" -rf "$SITES_DIR/$PREFIX.regions" \
-GL 1 -doGlf 2 \
-minMapQ 30 -minQ 20 \
-bam "$BAM_LIST" \
-remove_bads 1 \
-P 2

# settings:
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 3 takes major and minor allele from -sites file (for consistency)
# -bam list of bams to include
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -P K splits the analysis job over K nodes (but does not distribute I/O)
