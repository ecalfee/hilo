#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/
#SBATCH -J remapV4_RIMMA0625
#SBATCH -o /home/ecalfee/hilo/slurm-log/remapV4_RIMMA0625_%A_%a.out
#SBATCH -t 5:00
#SBATCH --mem=2G

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# this script takes a bam with reads aligned to v3 -> fastq -> maps with bwa to v4 genome
BAM_IN="landraces_fromLi/original/RIMMA0625_Andean.IndelRealigned"
ID="RIMMA0625"
DIR_OUT="landraces_fromLi/original"

# Submit part 1 slurm job: bam to fastq
echo "submitting script for BAM -> fastq for "${ID}
slurm1_message=$(sbatch --export="BAM_IN_PREFIX=${BAM_IN},FASTQ_OUT_PREFIX=${BAM_IN},ALL" ../scripts/bam2fastq.sh)
echo ${slurm1_message}

# Extract job ID from first slurm submission
if ! echo ${slurm1_message} | grep -q "[1-9][0-9]*$"; then
   echo "Job(s) submission failed."
   echo ${slurm1_message}
   exit 1
else
   job1_id=$(echo ${slurm1_message} | grep -oh "[1-9][0-9]*$")
fi

# Submit 2nd script which will launch only if job1 is successfully completed:
echo "Now submitting script to run bwa mem for "${ID}
sbatch --depend=afterok:${job1_id} \
-t 6-00:00:00 \
--mem=40G \
--export="ID=${ID},FASTQ1=${BAM_IN}_1.fq,FASTQ2=${BAM_IN}_2.fq,DIR_OUT=${DIR_OUT},ALL" \
../scripts/map2maizeAPGv4.sh
# above, I increased from default time and memory for bwa due to large file size (high coverage)

echo "all done!"
