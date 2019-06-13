#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/within_ancestry_diversity
#SBATCH -J freqLocus
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcAlleleFreq_onelocus_%A_%a.out
#SBATCH -t 12:00:00
#SBATCH --mem=8G
#SBATCH -n 1

# to run, e.g.
# sbatch --array=0-7 --export="PREFIX=pass2_alloMAIZE_parviglumis_inv4m_allSNPs,DIR_POPS=results/inv4m/pass2_alloMAIZE" calcAlleleFreq_onelocus.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset


REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
i=$SLURM_ARRAY_TASK_ID
pops=($(cat ${DIR_POPS}/pops.list))
POP=${pops["$i"]}
DIR_OUT="results/allele_freq/${PREFIX}/${POP}"
DIR_SCRATCH="/scratch/ecalfee/${PREFIX}/${POP}_alleleFreqCounts"


# load angsd v9.21 from bio
module load bio

# make directory to store output (if doesn't yet exist)
mkdir -p "${DIR_OUT}"
mkdir -p "${DIR_SCRATCH}"

echo "calculating allele freq. at variant sites region " $REGION_I "for pop "$POP

angsd -out "${DIR_SCRATCH}/whole_genome" \
-ref "$REF" \
-rf "${DIR_POPS}/regions.txt" \
-bam "${DIR_POPS}/${POP}.bams" \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 3 \
-sites "../global_ancestry/results/thinnedSNPs/${PREFIX}/whole_genome.var.sites" \
-doCounts 1 \
-doMaf 8 \
-P 1

echo "all done calculating pop allele frequencies! Now transfering files from local to home directory"
rsync -avh --remove-source-files ${DIR_SCRATCH}/ ${DIR_OUT}/ # copies all contents of output directory over to the appropriate home directory & cleans up scratch dir.
echo "results copied to home output directory: "${DIR_OUT}




# settings:
# -rf specifies regions file
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 3: takes major & minor allele from sites file
# -sites var.sites file has 4 tab separated columns: chrom pos major minor
# -bam list of bams to include
# -doMaf 8 is an unbiased estimator of pop freq using straight counts and no HWE assumption; -doMaf 8 with -doCounts 1
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
