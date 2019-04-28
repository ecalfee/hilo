#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/local_ancestry
#SBATCH -J getVCF4Allo
#SBATCH -o /home/ecalfee/hilo/slurm-log/getVCF4Allo_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=32G
#SBATCH --array=1-10

# to run: sbatch --export=PREFIX=pas2_alloMAIZE getVCFforAllo.sh

# note: number of threads should not exceed memory requested divided by 8G/thread

# array is the chromosome #
POP="allo.maize"
i=$SLURM_ARRAY_TASK_ID
DIR_SITES="results/thinnedSNPs/$PREFIX"
DIR_OUT="results/counts/$PREFIX"
DIR_SCRATCH="/scratch/ecalfee/vcfAllo_chr"${i}
BAMS="../samples/alloMAIZE_bams.list"
n_threads=2

# minimum and maximum individual depth filters. max will be set in export variables
MAX_DEPTH_IND=62 # 2x the mean for the highest coverage individual
MIN_DEPTH_IND=4

# this script takes in a sites file and calls genotypes for allopatric individuals
# from high coverage maize, outputing a VCF

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load angsd latest v9.21 using bio module
module load bio

# make temporary scratch directories
mkdir -p $DIR_SCRATCH
mkdir -p $DIR_OUT

echo "calling genotypes and outputting a vcf chr"${i}" pop "${POP}
# (1) make VCF in ANGSD
angsd -out ${DIR_SCRATCH}/${POP}_chr${i} \
-r ${i} \
-sites ${DIR_SITES}/chr${i}.var.sites \
-bam $BAMS \
-remove_bads 1 \
-minMapQ 30 \
-minQ 20 \
-doMajorMinor 3 \
-doMaf 1 \
-GL 1 \
-doGeno -4 \
-doPost 1 \
-postCutoff 0.9 \
-doVCF 1 \
-doCounts 1 \
-setMaxDepthInd ${MAX_DEPTH_IND} \
-geno_minDepth ${MIN_DEPTH_IND} \
-P ${n_threads}

# (maybe should be more lenient (.9) posterior prob. for calling genotypes or later additionally sample 1 read for low coverage positions in these allopatric maize individuals .. I'm not sure_

# settings:
# -doCounts allows us to filter based on maximum and minimum read depths
# -r specifies which region to work on
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
# -doVCF 1  outputs a vcf
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)

# copy results back over
echo "all done making allopatric maize VCF! Now transfering files from local to home directory"
rsync -avh --remove-source-files ${DIR_SCRATCH}/ ${DIR_OUT}/ # copies all contents of output directory over to the appropriate home directory & cleans up scratch dir.
echo "results copied to home output directory: "${DIR_OUT}

echo "all done!"
