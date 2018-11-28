#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J geneOverlap
#SBATCH -o /home/ecalfee/hilo/slurm-log/geneOverlapWindowsCM_%A_%a.out
#SBATCH -t 10:00:00
#SBATCH --mem=8G


DIR_WINDOWS="geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/windows0.1cM"
CDS_FILE="refMaize/geneAnnotations/CDS_merged.sh"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load software needed
module load bio

# first concatenate windows from all chromosomes, in order, into one .bed file: whole_genome.bed
echo "concatenating diff. chromosomes cM windows around SNPs to whole_genome.bed"
cat $(ls ${DIR_WINDOWS}/chr*.bed | sort --version-sort) > ${DIR_WINDOWS}/whole_genome.bed

# then calculate overlap:
echo "calculating gene overlap"
bedtools coverage -sorted -a ${DIR_WINDOWS}/whole_genome.bed -b ${CDS_FILE} > ${DIR_WINDOWS}/gene_overlap.bed
echo 'all done!'
