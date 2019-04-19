#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/global_ancestry
#SBATCH -J GL4var
#SBATCH -o /home/ecalfee/hilo/slurm-log/GL4variantsites_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=24G
#SBATCH -n 1
#SBATCH --array=1-10

# if you run on bigmemh or high2, need to limit array, e.g. --array=0-425%8
# %k ensures only k jobs max run at one time

# to run: global_ancestry$ sbatch --export=PREFIX=duplicates,SNPs=hilo_alloMAIZE_MAIZE4LOW_RIMMA0625_small GL_for_variant_sites.sh

# some VARIABLES
DIR_OUT=results/thinnedSNPs/"$PREFIX"
BAM_IN=../samples/"$PREFIX"_bams.list
i="$SLURM_ARRAY_TASK_ID"
REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
DIR_SCRATCH="/scratch/ecalfee/GL4_var_sites_chr_$i"
DIR_SITES=results/thinnedSNPs/"$SNPs"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# angsd v9.21 is in module bio
module load bio

# make directory to store output (if doesn't yet exist)
mkdir -p $DIR_OUT
mkdir -p $DIR_SCRATCH

# apply filtering with SAMtools & PICARD
echo "calculating GL's for variant sites chr "$i "for BAMS in "$BAM_IN
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) For each chromosomal region individually, find variant sites
angsd -out $DIR_SCRATCH/GL_chr$i \
-r $i \
-ref $REF \
-bam $BAM_IN \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 3 \
-sites ${DIR_SITES}/chr${i}.var.sites \
-GL 1 -doGlf 2 \
-P 1

echo "all done identifying variant sites! Now transfering files from local to home directory"
rsync -avh --remove-source-files ${DIR_SCRATCH}/ ${DIR_OUT}/ # copies all contents of output directory over to the appropriate home directory & cleans up scratch dir.
echo "results copied to home output directory: "${DIR_OUT}


# settings:
# -r specifies a region to work on (here: a chromosome)
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 3: use sites file to get major/minor allele
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
