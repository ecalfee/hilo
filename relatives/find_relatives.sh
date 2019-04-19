#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/relatives
#SBATCH -J relate
#SBATCH -o /home/ecalfee/hilo/slurm-log/findRelatives_%A_%a.out
#SBATCH -t 24:00:00
#SBATCH --mem=48G
#SBATCH -n 1
#SBATCH --export="PREFIX=specify_bam_set_here,REGIONS=regions_file_here,ALL"


# to run
# sbatch --export="PREFIX=duplicates,REGIONS=N1000.L100.regions,N_IND=74,ALL" find_relatives.sh


# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

DIR_OUT=results/"$REGIONS"
REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
BAM_IN=../samples/"$PREFIX"_bams.list

# load angsd v9.21 from bio
module load bio

# make directory to store output (if doesn't yet exist)
mkdir -p ${DIR_OUT}

echo "first calc likelihood of all 10 genotypes for each individual"

angsd -out ${DIR_OUT}/${PREFIX} \
-rf ../data/refMaize/random_regions/${REGIONS}.regions \
-ref "$REF" \
-bam "$BAM_IN" \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doGlf 1 \
-GL 1 \
-P 1

echo "all done calculating genotype likelihoods; now getting genotype fractions for each individual"
misc/ibs -f ${DIR_OUT}/${PREFIX}.glf.gz -nInd $N_IND -o all
# The output file is all.ibs

echo "now estimating the 10x10 genotype fraction matrix for all pairs (very slow)"
misc/ibs -f ${DIR_OUT}/${PREFIX}.glf.gz -nInd $N_IND -allpairs 1 -o all
# The output file is all.ibspair

# settings:
# -rf specifies regions file
# -remove_bads removes reads with flags like duplicates
# -bam list of bams to include
# -GL 1: use samtools genotype likelihood method
# -doGlf 1: prints all 10 possible genotypes and likelihoods
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
