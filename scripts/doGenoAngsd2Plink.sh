#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J doGeno2Plink
#SBATCH -o /home/ecalfee/hilo/slurm-log/doGeno2Plink_%j_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH --array=0-425
#SBATCH --export=POP=maize.allo.4Low16,DIR_POPS=pass1_bam_pops,DIR_REGIONS=refMaize/divide_5Mb,DIR_SITES=var_sites/merged_pass1_all_alloMaize4Low_16,DIR_OUT=geno/merged_pass1_all_alloMaize4Low_16

# to run
# sbatch --export=POP=maize.allo.4Low16,DIR_POPS=etc... doGenoAngsd2Plink.sh

# %k ensures only k jobs max run at one time, e.g. --array=0-425%8 runs 8 at a time

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

REGION_I=$SLURM_ARRAY_TASK_ID

# load angsd -- don't load -- updated version 9.20 is local
#module load angsd

# make directory to store output (if doesn't yet exist)
mkdir -p ${DIR_OUT}/${POP}

echo "inferring genotypes at variant sites region " $REGION_I "for pop "$POP

angsd -out ${DIR_OUT}/${POP}/region_${REGION_I} \
-rf ${DIR_REGIONS}/region_${REGION_I}.txt \
-ref refMaize/AGPv4.fa \
-bam ${DIR_POPS}/${POP}.list \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 3 \
-sites ${DIR_SITES}/region_${REGION_I}.var.sites \
-doMaf 1 \
-GL 1 \
-doGeno -4 \
-doPost 1 \
-postCutoff 0.95 \
-doPlink 2 \
-P 1

echo "all done!"

# settings:
# -r specifies which region to work on; -rf is the regions file
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 3: takes major & minor allele from sites file
# -sites var.sites file has 4 tab separated columns: chrom pos major minor
# -bam list of bams to include
# -GL 1: use samtools genotype likelihood method
# -doMaf 1 use fixed minor and major alleles. assumes HWE. alternative would be to use straight counts and no HWE assumption; -doMaf 8 with -doCounts 1
# (MAF is used as the prior for the genotype posterior)
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# -doGeno -4: suppresses geno output (see ANGSD documentation plink page)
# -doPost 1: calculate genotype posterior using MAF as prior
# -postCutoff sets the posterior confidence below which no genotype is called (set instead as missing)
# I do not set a separate cutoff for # of reads needed to call a genotype -geno_minDepth 4 (requires -doCounts)
# -doPlink 2 outputs genotypes in a plink formatted .tfam .tped files
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)