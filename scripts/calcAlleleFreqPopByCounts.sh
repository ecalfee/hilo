#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J countFreq
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcAlleleFreqPopByCounts_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH --array=0-425
#SBATCH --export="POP=specify_pop_here,ALL"

DIR_POPS="pass1_bam_pops"
DIR_SITES="geno_lik/merged_pass1_all_alloMaize4Low_16/allVar_depthFilt"
DIR_POP_FREQS="geno_lik/merged_pass1_all_alloMaize4Low_16/allVar_depthFilt/popFreqs"

# to run
# sbatch --export="POP=pop363,ALL" calcAlleleFreqPopByCounts.sh
# sbatch --export="POP=maize.allo.4Low16,ALL" calcAlleleFreqPopByCounts.sh

# %k ensures only k jobs max run at one time, e.g. --array=0-425%8 runs 8 at a time

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset
REGION_I=$SLURM_ARRAY_TASK_ID
DIR_SCRATCH=/scratch/ecalfee/${POP}_alleleFreqCounts/region_${REGION_I}
DIR_OUT=${DIR_POP_FREQS}/${POP}

# load angsd v9.21 from bio
module load bio

# make directory to store output (if doesn't yet exist)
mkdir -p ${DIR_OUT}
mkdir -p ${DIR_SCRATCH}

echo "calculating allele freq. at variant sites region " $REGION_I "for pop "$POP

angsd -out ${DIR_SCRATCH}/region_${REGION_I} \
-ref refMaize/AGPv4.fa \
-bam ${DIR_POPS}/${POP}.list \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 3 \
-sites ${DIR_SITES}/region_${REGION_I}.var.sites \
-doCounts 1 \
-doMaf 8 \
-P 1

echo "all done calculating pop allele frequencies! Now transfering files from local to home directory"
rsync -avh --remove-source-files ${DIR_SCRATCH}/ ${DIR_OUT}/ # copies all contents of output directory over to the appropriate home directory & cleans up scratch dir.
echo "results copied to home output directory: "${DIR_OUT}




# settings:
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 3: takes major & minor allele from sites file
# -sites var.sites file has 4 tab separated columns: chrom pos major minor
# -bam list of bams to include
# -doMaf 8 is an unbiased estimator of pop freq using straight counts and no HWE assumption; -doMaf 8 with -doCounts 1
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
