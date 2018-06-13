#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J varAngsdSymp
#SBATCH -o /home/ecalfee/hilo/slurm-log/varAngsdSymp_%j_%A_%a.out
#SBATCH -t 10:00:00
#SBATCH --mem=15G
#SBATCH -n 4

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load angsd
module load angsd

# make directory to store output (if doesn't yet exist)
mkdir -p var_sites/pass1/sympatric/

# apply filtering with SAMtools & PICARD
echo "calling variants using ANGSD on sympatric BAMS for hilo CHR$SLURM_ARRAY_TASK_ID"
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) For each chromosome individually, find variant sites
angsd -out var_sites/pass1/sympatric/region_$SLURM_ARRAY_TASK_ID \
-r $(cat refMaize/divide_50Mb/region_$SLURM_ARRAY_TASK_ID.txt) \
-rf refMaize/divide_50Mb/region_$SLURM_ARRAY_TASK_ID.txt
-doMajorMinor 2 \
-bam pass1_bam.onlySympatric.list \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doCounts 1 -minMaf 0.05 -doMaf 8 \
-minInd 50 \
-P 4

# settings:
# -r specifies which region to work on
# -rf specifies the regions in a file
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 4: pre-specify major allele from reference genome and infer minor allele from genotype likelihood
# -bam list of bams to include (all newly sequenced sympatric mexicana & maize pops)
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -doMaf 2 : output minor allele freq
# an alternative would be to use AlleleCounts method (-doMaf 8) with -doCounts
# which doesn't consider base quality and ignores all reads with non-target alleles
# but also doesn't rely on HW equlibrium assumptions to get ML allele freq from genotype freqs
# -minMaf x: and then do a cutoff to only include variant sites with >x diff. from reference allele
# -minInd N: only keep sites with information (at least one read) from N individuals
# -P 4 means use 4 threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
