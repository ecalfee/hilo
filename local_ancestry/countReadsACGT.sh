#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/local_ancestry
#SBATCH -J countACGT
#SBATCH -o /home/ecalfee/hilo/slurm-log/countACGT_%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --mem=8G

# to run: local_ancestry$ sbatch --array=0-217 --export="CHR=1,PREFIX=pass2_alloMAIZE,ALL" countReadsACGT.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# this script takes in a sites file and outputs an angsd counts file with
# columns totA totC totG totT as well as a .pos.gz file with all positions with data
# for a low-coverage hilo individual

MAX_DEPTH_IND=8
REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
i=$SLURM_ARRAY_TASK_ID # which hilo sample from list file
ids=($(cat ../samples/"$PREFIX"_IDs.list))
ID=${ids["$i"]} # sample ID
bams=($(cat ../samples/"$PREFIX"_bams.list))
BAM=${bams["$i"]} # path to individual bam

DIR_SITES="results/thinnedSNPs/$PREFIX/"
DIR_OUT="results/countsACGT/$PREFIX/$ID"
DIR_SCRATCH="/scratch/ecalfee/countACGT_"$ID"_chr$i"

# load modules
# includes angsd v9.21
module load bio

# make output and scratch directories
mkdir -p "$DIR_OUT"
mkdir -p "$DIR_SCRATCH"

# make file with single focus bam listed
mkdir -p "../samples/"$PREFIX"_list_individuals"
echo "$BAM" > "../samples/"$PREFIX"_list_individuals/"$ID"_bams.list"

echo "counting reads supporting ACGT alleles for chr"${CHR}" and individual "${ID}
# counts reads for sites in sites file that pass basic quality filtering
angsd -out "$DIR_SCRATCH/chr$CHR" \
-ref "$REF" \
-bam "../samples/"$PREFIX"_list_individuals/"$ID"_bams.list" \
-minQ 20 -minMapQ 30 \
-doCounts 1 -dumpCounts 3 \
-setMaxDepthInd "$MAX_DEPTH_IND" \
-remove_bads 1 \
-r "$CHR" \
-sites "$DIR_SITES/chr$CHR.var.sites"

# transfer output files to home directory
echo "all done counting ACGT! Now transfering files from local to home directory"
rsync -avh --remove-source-files "$DIR_SCRATCH/" "$DIR_OUT/" # copies all contents of output directory over to the appropriate home directory & cleans up scratch dir.
echo "results copied to home output directory: "${DIR_OUT}

echo "all done!"
