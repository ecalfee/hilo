#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/within_ancestry_diversity
#SBATCH -J piFst
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcPiFstAngsd_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH -n 2
#SBATCH --mem=4G

# this script calculate pi and fst for outlier regions from a file
# to run:
# within_ancestry_diversity$ sbatch --array=0-25 --export=OUTLIERS=../ZAnc/results/pass1_maize/outliers.bed,PREFIX=mex2,SUBDIR_OUT=piFst_outliers piFstAngsd_outliers.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# loads angsd
module load bio

i=$SLURM_ARRAY_TASK_ID # 0 indexed!


LIST_CHR=($(cut -f1 $OUTLIERS))
CHR="${LIST_CHR["$i"]}"
LIST_START=($(cut -f2 $OUTLIERS))
START="${LIST_START["$i"]}"
LIST_END=($(cut -f3 $OUTLIERS))
END="${LIST_END["$i"]}"
POPS=(symp.mexicana symp.maize allo.mexicana pop18 pop19 pop20 pop21 pop22 pop23 pop24 pop25 pop26 pop27 pop28 pop29 pop30 pop31 pop33 pop34 pop35 pop360 pop361 pop362 pop363 pop365 pop366 pop367 pop368 pop369 pop370 pop371 pop372 pop373 pop374)

DIR_OUT="results/$SUBDIR_OUT/$PREFIX/chr${CHR}_${START}-${END}"
REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"

echo "outlier region: $CHR:$START-$END"

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
-bam "results/input/$PREFIX/pops/$POP.list" \
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
-bam "results/input/$PREFIX/pops/$POP.list" \
-remove_bads 1 -minMapQ 30 -minQ 20 \
-GL 1 \
-P 2

echo "summarizing thetas for region and windows"
thetaStat do_stat "$DIR_OUT"/"$POP".thetas.idx -outnames "$DIR_OUT"/"$POP".thetasAll
thetaStat do_stat "$DIR_OUT"/"$POP".thetas.idx -win 5000 -step 1000 -outnames "$DIR_OUT"/"$POP".thetasWindows

done

# calculate pairwise 2D SFS
echo "calculating pairwise 2D SFS for each pair of pops"
realSFS "$DIR_OUT/symp.maize.saf.idx" "$DIR_OUT/symp.mexicana.saf.idx" -P 2 > "$DIR_OUT/symp.maize-symp.mexicana.sfs"
realSFS "$DIR_OUT/symp.mexicana.saf.idx" "$DIR_OUT/allo.mexicana.saf.idx" -P 2 > "$DIR_OUT/symp.mexicana-allo.mexicana.sfs"
realSFS "$DIR_OUT/symp.maize.saf.idx" "$DIR_OUT/allo.mexicana.saf.idx" -P 2 > "$DIR_OUT/symp.maize-allo.mexicana.sfs"

echo "now calculating FST for each pair"
# calculate FST (-whichFst 1 indicates Hudson/Bhatia2013 estimator for Fst)
# sympatric maize and mexicana
realSFS fst index "$DIR_OUT/symp.maize.saf.idx" "$DIR_OUT/symp.mexicana.saf.idx" \
-sfs "$DIR_OUT/symp.maize-symp.mexicana.sfs" \
-whichFst 1 \
-fstout "$DIR_OUT/symp.maize-symp.mexicana"
# future analyses should rely on 'weighted' Fst (i.e. the ratio of expectations across loci).
realSFS fst stats "$DIR_OUT/symp.maize-symp.mexicana.fst.idx" > "$DIR_OUT/symp.maize-symp.mexicana.fst.stats"
# summary fst
realSFS fst stats2 "$DIR_OUT/symp.maize-symp.mexicana.fst.idx" > "$DIR_OUT/symp.maize-symp.mexicana.fstAll"
# windows fst
realSFS fst stats2 "$DIR_OUT/symp.maize-symp.mexicana.fst.idx" -win 5000 -step 1000 > "$DIR_OUT/symp.maize-symp.mexicana.fstWindows.stats"

# sympatric and allopatric mexicana
realSFS fst index "$DIR_OUT/symp.mexicana.saf.idx" "$DIR_OUT/allo.mexicana.saf.idx" \
-sfs "$DIR_OUT/symp.mexicana-allo.mexicana.sfs" \
-whichFst 1 \
-fstout "$DIR_OUT/symp.mexicana-allo.mexicana"
# future analyses should rely on 'weighted' Fst (i.e. the ratio of expectations across loci).
realSFS fst stats "$DIR_OUT/symp.mexicana-allo.mexicana.fst.idx" > "$DIR_OUT/symp.mexicana-allo.mexicana.fst.stats"
# summary fst
realSFS fst stats2 "$DIR_OUT/symp.mexicana-allo.mexicana.fst.idx" > "$DIR_OUT/symp.mexicana-allo.mexicana.fstAll"
# windows fst
realSFS fst stats2 "$DIR_OUT/symp.mexicana-allo.mexicana.fst.idx" -win 5000 -step 1000 > "$DIR_OUT/symp.mexicana-allo.mexicana.fstWindows.stats"

# sympatric maize and allopatric mexicana
realSFS fst index "$DIR_OUT/symp.maize.saf.idx" "$DIR_OUT/allo.mexicana.saf.idx" \
-sfs "$DIR_OUT/symp.maize-allo.mexicana.sfs" \
-whichFst 1 \
-fstout "$DIR_OUT/symp.maize-allo.mexicana"
# future analyses should rely on 'weighted' Fst (i.e. the ratio of expectations across loci).
realSFS fst stats "$DIR_OUT/symp.maize-allo.mexicana.fst.idx" > "$DIR_OUT/symp.maize-allo.mexicana.fst.stats"
# summary fst
realSFS fst stats2 "$DIR_OUT/symp.maize-allo.mexicana.fst.idx" > "$DIR_OUT/symp.maize-allo.mexicana.fstAll"
# windows fst
realSFS fst stats2 "$DIR_OUT/symp.maize-allo.mexicana.fst.idx" -win 5000 -step 1000 > "$DIR_OUT/symp.maize-allo.mexicana.fstWindows.stats"

echo "all done!"

# options
# basic quality filtering for reads
# anc polarizes SFS by reference genome
# underFlowProtect is necessary for large #s of bams
