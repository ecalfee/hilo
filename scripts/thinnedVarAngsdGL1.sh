#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J thinGL1Angsd
#SBATCH -o /home/ecalfee/hilo/slurm-log/thinnedGL1Angsd_%A_%a.out
#SBATCH -t 24:00:00
#SBATCH --mem=32G
#SBATCH -n 1
#SBATCH --array=1-10

# some VARIABLES
i=$SLURM_ARRAY_TASK_ID #chromosome
DIR_SITES="geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedPCA"
DIR_OUT=$DIR_SITES
BAM_IN="merged_bam.pass1_all.alloMaize4Low_16.list"
DIR_SCRATCH="/scratch/ecalfee/thinnedPCA_chr"${i}
MAX_DEPTH_IND=63

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
echo "calling variants and GL using ANGSD on BAMS for hilo genomic regions "$SLURM_ARRAY_TASK_ID
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) For each chromosomal region individually, find variant sites
angsd -out $DIR_SCRATCH/chr_$i \
-r ${i} \
-sites ${DIR_SITES}/chr${i}.var.sites \
-bam $BAM_IN \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 3 \
-GL 1 -doGlf 2 \
-P 1 \
-setMaxDepthInd ${MAX_DEPTH_IND}

# took out ref -ref refMaize/AGPv4.fa \

echo "all done identifying variant sites! Now transfering files from local to home directory"
rsync -avh --remove-source-files ${DIR_SCRATCH}/ ${DIR_OUT}/ # copies all contents of output directory over to the appropriate home directory & cleans up scratch dir.
echo "results copied to home output directory: "${DIR_OUT}


# settings:
# -r specifies chromosome to work on
# -setMaxDepthInd excludes individuals at any site where they exceed ind max depth
# -sites var.sites file has 4 tab separated columns: chrom pos major minor
# -remove_bads removes reads with flags like duplicates
# -bam list of bams to include (all newly sequenced allopatric mex. and sympatric mexicana & maize pops)
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -doMajorMinor 3: takes major & minor allele from sites file
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
