#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/variant_sites
#SBATCH -J idSNPs
#SBATCH -o /home/ecalfee/hilo/slurm-log/idSNPs_%A_%a.out
#SBATCH -t 10:00:00
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH --array=0-425

# if you run on bigmemh or high2, need to limit array, e.g. --array=0-425%8
# %k ensures only k jobs max run at one time

# to run: variant_sites$ sbatch --export=PREFIX=hilo_alloMAIZE_MAIZE4LOW,MAX_DEPTH=880 id_SNPs.sh

# some VARIABLES
DIR_OUT=results/"$PREFIX"
BAM_IN=../samples/"$PREFIX"_bams.list
i="$SLURM_ARRAY_TASK_ID"
DIR_REGIONS="../data/refMaize/divide_5Mb"
REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
DIR_SCRATCH="/scratch/ecalfee/id_SNPs_region_$i"

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
echo "calling variants and GL using ANGSD on BAMS for hilo genomic region "$i
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) For each chromosomal region individually, find variant sites
angsd -out $DIR_SCRATCH/region_$i \
-rf $DIR_REGIONS/region_$i.txt \
-ref $REF \
-bam $BAM_IN \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 2 \
-doCounts 1 -minMaf 0.05 -doMaf 8 \
-GL 1 -doGlf 2 \
-minInd 50 \
-P 1 \
-setMaxDepth ${MAX_DEPTH}

echo "all done identifying variant sites! Now transfering files from local to home directory"
rsync -avh --remove-source-files ${DIR_SCRATCH}/ ${DIR_OUT}/ # copies all contents of output directory over to the appropriate home directory & cleans up scratch dir.
echo "results copied to home output directory: "${DIR_OUT}


# settings:
# -rf specifies in a file which region to work on
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 2: infer major and minor from allele counts
# -bam list of bams to include (all newly sequenced allopatric mex. and sympatric mexicana & maize pops)
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# we use AlleleCounts method for MAF (-doMaf 8) with -doCounts
# which doesn't consider base quality and ignores all reads with non-target alleles
# but also doesn't rely on HW equlibrium assumptions to get ML allele freq from genotype freqs
# -minMaf x: and then do a cutoff to only include variant sites with >x minor allele freq.
# -minInd N: only keep sites with information (at least one read) from N individuals
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
# -setMaxDepth filters out sites where total depth exceeds some threshold
