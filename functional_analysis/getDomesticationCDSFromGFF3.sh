#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J domestCDS
#SBATCH -o /home/ecalfee/hilo/slurm-log/getDomesticationCDSfromGFF_%A_%a.out
#SBATCH -t 1:00:00
#SBATCH --mem=8G

DIR="refMaize/geneAnnotations"
CDS_GFF=$DIR"/CDS_autosome_only.gff" # output of mergeCDSFromGFF3.sh
CHROM_ORDER_FILE="refMaize/Zea_mays.AFPv4.dna.chr.autosome.lengths"
DOMESTICATION_LIST="domestication/gene_model_translation_to_APGv4.txt"
CDS_DOMESTICATION=$DIR"/CDS_domestication_Hufford2012_merged.bed"
DIR_WINDOWS_CM="geno_lik/merged_pass1_all_alloMaize4Low_16/thinnedHMM/windows0.1cM"
DIR_WINDOWS_BP="refMaize/windows_10kb"

# note: domestication list of genes was acquired from:
# gene names from Hufford 2012's original paper to the names used for the maize v4 genome assembly:
# https://www.maizegdb.org/gene_center/gene# "Translate Gene Model IDs" to Zm00001d.2.
# Saved output file as gene_model_translation_to_APGv4.txt

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load programs
module load bio

# find domestication genes and gff3 overlap:
for i in $(cut -f 3 $DOMESTICATION_LIST | uniq | awk '$0 != "" && $0 != "Zm00001d.2" {print $0}'); \
do grep $i $CDS_GFF; done > $DIR/CDS_domestication_Hufford2012_only.gff

# sort by coordinates and merge overlapping or contiguous CDS
echo "sorting by coordinates and merging overlapping or contiguous CDS"
bedtools sort -g ${CHROM_ORDER_FILE} -i $DIR/CDS_domestication_Hufford2012_only.gff | \
bedtools merge  > $CDS_DOMESTICATION

# then calculate overlap:
echo "calculating domestication gene overlap with fixed cM windows around SNPs"
bedtools coverage -sorted -a ${DIR_WINDOWS_CM}/whole_genome.bed -b ${CDS_DOMESTICATION} \
> ${DIR_WINDOWS_CM}/domestication_Hufford2012_overlap.bed

echo "calculating domestication gene overlap with fixed bp windows genomewide"
bedtools coverage -sorted -a ${DIR_WINDOWS_BP}/whole_genome.bed -b ${CDS_DOMESTICATION} \
> ${DIR_WINDOWS_BP}/domestication_Hufford2012_overlap.bed
echo 'all done!'
