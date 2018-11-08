#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J calcPT
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcAlleleFreqPopsTransfer_%A_%a.out
#SBATCH -t 2:00:00
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH --array=0-425
#SBATCH --export="POP=specify_pop_here,ALL"

DIR_POPS="pass1_bam_pops"
DIR_REGIONS="refMaize/divide_5Mb"
DIR_SITES="geno_lik/merged_pass1_all_alloMaize4Low_16/allVar_depthFilt"

# to run
# sbatch --export="POP=pop363,ALL" calcAlleleFreqPop.sh
# sbatch --export="POP=maize.allo.4Low16,ALL" calcAlleleFreq.sh

# %k ensures only k jobs max run at one time, e.g. --array=0-425%8 runs 8 at a time

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset
REGION_I=$SLURM_ARRAY_TASK_ID
DIR_SCRATCH=/scratch/ecalfee/${POP}/region_${REGION_I}
DIR_OUT=${DIR_SITES}/${POP}

# load angsd v9.21 from bio
module load bio

# make directory to store output (if doesn't yet exist)
mkdir -p ${DIR_OUT}
mkdir -p ${DIR_SCRATCH}

echo "calculating allele freq. at variant sites region " $REGION_I "for pop "$POP

angsd -out ${DIR_SCRATCH}/region_${REGION_I} \
-rf ${DIR_REGIONS}/region_${REGION_I}.txt \
-ref refMaize/AGPv4.fa \
-bam ${DIR_POPS}/${POP}.list \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 3 \
-sites ${DIR_SITES}/region_${REGION_I}.var.sites \
-doMaf 1 \
-GL 1 \
-P 1

echo "all done calculating pop allele frequencies! Now transfering files from local to home directory"
rsync -avh ${DIR_SCRATCH}/ ${DIR_OUT}/ # copies all contents of output directory over to the appropriate home directory
echo "results copied to home output directory: "${DIR_OUT}




# settings:
# -r specifies which region to work on; -rf is the regions file
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 3: takes major & minor allele from sites file
# -sites var.sites file has 4 tab separated columns: chrom pos major minor
# -bam list of bams to include
# -GL 1: use samtools genotype likelihood method
# -doMaf 1 use fixed minor and major alleles. assumes HWE. alternative would be to use straight counts and no HWE assumption; -doMaf 8 with -doCounts 1
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
