#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/within_ancestry_diversity
#SBATCH -J piFst1
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcPiFst1locus_%A_%a.out
#SBATCH -t 24:00:00
#SBATCH -n 2
#SBATCH --mem=16G

# this script calculate pi and fst for one locus of interest
# using allopatric maize and allopatric mexicana homozygous for ancestry haplotypes
# to run:
# within_ancestry_diversity$ sbatch --export=PREFIX=pass2_alloMAIZE,LOCUS=inv4m,CHR=4,START=170771502,END=186951149 piFstAngsd_onelocus.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# loads angsd
module load bio
POPS=(allopatric_mexicana.mexicana_hap allopatric_maize.maize_hap)

DIR_OUT="results/$LOCUS/$PREFIX"
REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"

echo "outlier region $LOCUS: $CHR:$START-$END"

# make directory to store output (if doesn't yet exist)
mkdir -p "$DIR_OUT"

echo "finding site allele frequencies per pop"
for POP in ${POPS[@]}
do echo $POP
# steps:
# (0) Start with filtered BAM files and reference genome
# (1) Use all sites to estimate site allele frequency
angsd -out "$DIR_OUT/$POP" \
-anc "$REF" \
-fold 0 \
-underFlowProtect 1 \
-r "$CHR:$START-$END" \
-bam "results/$LOCUS/$PREFIX/$POP.bams" \
-remove_bads 1 -minMapQ 30 -minQ 20 \
-GL 1 \
-doSaf 1 \
-P 2

echo "done with SAF! calculating SFS"
realSFS "$DIR_OUT/$POP.saf.idx" -P 4 > "$DIR_OUT/$POP.sfs"

echo "done with SFS! calculating within-pop diversity 'thetas'"
angsd -out "$DIR_OUT/$POP" \
-anc "$REF" \
-doThetas 1 \
-doSaf 1 \
-pest "$DIR_OUT/$POP.sfs" \
-underFlowProtect 1 \
-r "$CHR:$START-$END" \
-bam "results/$LOCUS/$PREFIX/$POP.bams" \
-remove_bads 1 -minMapQ 30 -minQ 20 \
-GL 1 \
-P 2

echo "summarizing thetas for region and windows"
thetaStat do_stat "$DIR_OUT"/"$POP".thetas.idx -outnames "$DIR_OUT"/"$POP".thetasAll
thetaStat do_stat "$DIR_OUT"/"$POP".thetas.idx -win 5000 -step 1000 -outnames "$DIR_OUT"/"$POP".thetasWindows

done

# calculate pairwise 2D SFS
pop1="allopatric_mexicana.mexicana_hap"
pop2="allopatric_maize.maize_hap"


echo "calculating pairwise 2D SFS for each pair of pops"
realSFS "$DIR_OUT/$pop1.saf.idx" "$DIR_OUT/$pop2.saf.idx" -P 2 > "$DIR_OUT/$pop1-$pop2.sfs"

echo "now calculating FST for each pair"
# calculate FST (-whichFst 1 indicates Hudson/Bhatia2013 estimator for Fst)
# allopatric maize and mexicana
realSFS fst index "$DIR_OUT/$pop1.saf.idx" "$DIR_OUT/$pop2.saf.idx" \
-sfs "$DIR_OUT/$pop1-$pop2.sfs" \
-whichFst 1 \
-fstout "$DIR_OUT/$pop1-$pop2"
# future analyses should rely on 'weighted' Fst (i.e. the ratio of expectations across loci).
realSFS fst stats "$DIR_OUT/$pop1-$pop2.fst.idx" > "$DIR_OUT/$pop1-$pop2.fst.stats"
# summary fst
realSFS fst stats2 "$DIR_OUT/$pop1-$pop2.fst.idx" > "$DIR_OUT/$pop1-$pop2.fstAll"
# windows fst
realSFS fst stats2 "$DIR_OUT/$pop1-$pop2.fst.idx" -win 5000 -step 1000 > "$DIR_OUT/$pop1-$pop2.fstWindows.stats"

echo "all done!"

# options
# basic quality filtering for reads
# anc polarizes SFS by reference genome
# underFlowProtect is necessary for large #s of bams
