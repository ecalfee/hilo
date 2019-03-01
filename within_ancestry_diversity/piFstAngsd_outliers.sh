#!/bin/bash -l
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/within_ancestry_diversity
#SBATCH -J piFst
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcPiFstAngsd_%A_%a.out
#SBATCH -t 6:00:00
#SBATCH -n 4
#SBATCH --mem=8G

# this script calculate pi and fst for outlier regions from a file
# to run:
# within_ancestry_diversity$ sbatch --array=0-25 --export=OUTLIERS=../ZAnc/results/pass1_maize/outliers.bed,PREFIX=mex2 piFstAngsd_outliers.sh

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
END="${LIST_START["$i"]}"

DIR_OUT="results/piFst_outliers/$PREFIX/chr${CHR}_${START}-${END}"
REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"

echo "outlier region: $CHR:$START-$END"

# make directory to store output (if doesn't yet exist)
mkdir -p "$DIR_OUT"

echo "finding site allele frequencies per pop"
for POP in symp.mexicana symp.maize allo.mexicana
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
-P 4

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
-P 4

done

# calculate pairwise 2D SFS
echo "calculating pairwise 2D SFS for each pair of pops"
realSFS "$DIR_OUT/symp.mexicana.saf.idx" "$DIR_OUT/symp.maize.saf.idx" -P 4 > "$DIR_OUT/symp.maize-symp.mexicana.sfs"
realSFS "$DIR_OUT/symp.mexicana.saf.idx" "$DIR_OUT/allo.mexicana.saf.idx" -P 4 > "$DIR_OUT/symp.mexicana-allo.mexicana.sfs"
realSFS "$DIR_OUT/symp.maize.saf.idx" "$DIR_OUT/allo.mexicana.saf.idx" -P 4 > "$DIR_OUT/symp.maize-allo.mexicana.sfs"

echo "now calculating FST for each pair"
# calculate FST (-whichFst 1 indicates Hudson/Bhatia2013 estimator for Fst)
# sympatric maize and mexicana
realSFS fst index "$DIR_OUT/symp.mexicana.saf.idx" "$DIR_OUT/symp.maize.saf.idx" \
-sfs "$DIR_OUT/symp.maize-symp.mexicana.sfs" \
-whichFst 1 \
-fstout "$DIR_OUT/symp.maize-symp.mexicana"
# future analyses should rely on 'weighted' Fst (i.e. the ratio of expectations across loci).
realSFS fst stats "$DIR_OUT/symp.maize-symp.mexicana.fst.idx" > "$DIR_OUT/symp.maize-symp.mexicana.fst.stats"

# sympatric and allopatric mexicana
realSFS fst index "$DIR_OUT/symp.mexicana.saf.idx" "$DIR_OUT/allo.mexicana.saf.idx" \
-sfs "$DIR_OUT/symp.mexicana-allo.mexicana.sfs" \
-whichFst 1 \
-fstout "$DIR_OUT/symp.mexicana-allo.mexicana"
# future analyses should rely on 'weighted' Fst (i.e. the ratio of expectations across loci).
realSFS fst stats "$DIR_OUT/symp.mexicana-allo.mexicana.fst.idx" > "$DIR_OUT/symp.mexicana-allo.mexicana.fst.stats"

# sympatric maize and allopatric mexicana
realSFS fst index "$DIR_OUT/symp.maize.saf.idx" "$DIR_OUT/allo.mexicana.saf.idx" \
-sfs "$DIR_OUT/symp.maize-allo.mexicana.sfs" \
-whichFst 1 \
-fstout "$DIR_OUT/symp.maize-allo.mexicana"
# future analyses should rely on 'weighted' Fst (i.e. the ratio of expectations across loci).
realSFS fst stats "$DIR_OUT/symp.maize-allo.mexicana.fst.idx" > "$DIR_OUT/symp.maize-allo.mexicana.fst.stats"

echo "all done!"


# options
# basic quality filtering for reads
# anc polarizes SFS by reference genome
# underFlowProtect is necessary for large #s of bams
