#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J countReads
#SBATCH -o /home/ecalfee/hilo/slurm-log/countReadsRefAlt_%j_%A_%a.out
#SBATCH -t 8:00:00
#SBATCH --mem=8G
#SBATCH --array=1-200
#SBATCH --export=CHR=1

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# this script takes in a VCF sites file (from high coverage maize, produced by plink)
# and outputs a read counts file for the ref and alt allele for a low-coverage hilo individual

# array is the hilo ID #
id=$SLURM_ARRAY_TASK_ID
VCF_POP="maize.allo.4Low16" # note the VCF is from plink (not ANGSD directly -- can't read that format)
dir="var_sites/merged_pass1_all_alloMaize4Low_16/thinnedHMM/"

# load modules
# load GATK and dependencies
module load R
module load maven
module load java
module load GATK
# -h to see help
#java -Xmx6g -jar /share/apps/GATK-3.6/GenomeAnalysisTK.jar -T ASEReadCounter -h
echo "counting ref and alt reads per individual in bams for chr"${CHR}" hilo id = "${id}
java -Xmx7g -jar /share/apps/GATK-3.6/GenomeAnalysisTK.jar \
-T ASEReadCounter \
--sitesVCFFile ${dir}/${VCF_POP}_chr${CHR}.vcf \
--out ${dir}/hilo_${id}_chr${CHR}.csv \
--minBaseQuality 20 \
--minMappingQuality 30 \
-I filtered_bam/hilo_${id}.sort.dedup.baq.rg.bam \
-L $CHR

echo "all done!"
