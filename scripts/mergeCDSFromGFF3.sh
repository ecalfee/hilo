#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J mergeCDS
#SBATCH -o /home/ecalfee/hilo/slurm-log/mergeCDSfromGFF_%A_%a.out
#SBATCH -t 1:00:00
#SBATCH --mem=8G

DIR="refMaize/geneAnnotations"
ORIGINAL_GFF=$DIR"/Zea_mays.B73_RefGen_v4.41.chr.gff3"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load programs
module load bio

# limit intervals to only CDS
echo "filter to only CDS within original GFF3"
gffread -C -o- $ORIGINAL_GFF | awk '$3=="CDS" {print $0}' > $DIR/CDS_only.gff

# sort by coordinates and merge overlapping or contiguous CDS
echo "sorting by coordinates and merging overlapping or contiguous CDS"
bedtools sort -i $DIR/CDS_only.gff | bedtools merge  > $DIR/CDS_merged.bed

echo 'all done!'
