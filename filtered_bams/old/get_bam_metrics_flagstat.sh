#!/bin/bash
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/filtered_bams
#SBATCH -J flagStat
#SBATCH -o /home/ecalfee/hilo/slurm-log/getBamMetricsFlagstat_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=8G
#SBATCH --array=0

# to run: sbatch --array=0-296 --export=PREFIX=hilo_alloMAIZE_MAIZE4LOW get_bam_metrics_flagstat.sh

# START WITH ZERO FOR ARRAY INDEXES (!)

# this script runs samtools flagstat on a set of bam files

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
# load bio for samtools
module load bio

i=$SLURM_ARRAY_TASK_ID
ids=($(cat ../samples/${PREFIX}_IDs.list))
ID=${ids["$i"]}
bams=($(cat ../samples/${PREFIX}_bams.list))
BAM=${bams["$i"]}

# make output directory
mkdir -p metrics

echo "getting summary stats from samtools flagstat."
echo "ID: $ID"
echo "bam: $BAM"
samtools flagstat "$BAM" > "metrics/$ID.flagstat"

echo "all done!"
