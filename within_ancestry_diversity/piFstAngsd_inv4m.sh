#!/bin/bash -l
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/within_ancestry_diversity
#SBATCH -J inv4m_pi
#SBATCH -o /home/ecalfee/hilo/slurm-log/calcPiFstinv4m_%A_%a.out
#SBATCH -t 24:00:00
#SBATCH -n 2
#SBATCH --mem=16G

# this script calculate pi and fst for inv4m
# to run:
# within_ancestry_diversity$ sbatch --export=PREFIX=pass2_alloMAIZE piFstAngsd_inv4m.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# loads angsd
module load bio
#inv4m	4:171771502-185951149 but I add a 1Mb buffer on either side of the inversion
CHR=4
START=170771502
END=186951149
POPS=(allopatric_mexicana.mexicana_hap allopatric_maize.maize_hap parviglumis sympatric_maize.mexicana_hap sympatric_maize.maize_hap sympatric_mexicana.mexicana_hap sympatric_mexicana.maize_hap)

DIR_OUT="results/inv4m/$PREFIX"
REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"

echo "outlier region inv4m: $CHR:$START-$END"

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
-bam "results/inv4m/$PREFIX/$POP.bams" \
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
-bam "results/inv4m/$PREFIX/$POP.bams" \
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
pop3="parviglumis"

echo "calculating pairwise 2D SFS for each pair of pops"
realSFS "$DIR_OUT/$pop1.saf.idx" "$DIR_OUT/$pop2.saf.idx" -P 2 > "$DIR_OUT/$pop1-$pop2.sfs"
realSFS "$DIR_OUT/$pop3.saf.idx" "$DIR_OUT/$pop1.saf.idx" -P 2 > "$DIR_OUT/$pop3-$pop1.sfs"
realSFS "$DIR_OUT/$pop3.saf.idx" "$DIR_OUT/$pop2.saf.idx" -P 2 > "$DIR_OUT/$pop3-$pop1.sfs"

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

# parviglumis and allopatric maize
realSFS fst index "$DIR_OUT/$pop3.saf.idx" "$DIR_OUT/$pop1.idx" \
-sfs "$DIR_OUT/$pop3-$pop1.sfs" \
-whichFst 1 \
-fstout "$DIR_OUT/$pop3-$pop1"
# future analyses should rely on 'weighted' Fst (i.e. the ratio of expectations across loci).
realSFS fst stats "$DIR_OUT/$pop3-$pop1.fst.idx" > "$DIR_OUT/$pop3-$pop1.fst.stats"
# summary fst
realSFS fst stats2 "$DIR_OUT/$pop3-$pop1.fst.idx" > "$DIR_OUT/$pop3-$pop1.fstAll"
# windows fst
realSFS fst stats2 "$DIR_OUT/$pop3-$pop1.fst.idx" -win 5000 -step 1000 > "$DIR_OUT/$pop3-$pop1.fstWindows.stats"

# parviglumis and allopatric mexicana
realSFS fst index "$DIR_OUT/$pop3.saf.idx" "$DIR_OUT/$pop2.saf.idx" \
-sfs "$DIR_OUT/$pop3-$pop2.sfs" \
-whichFst 1 \
-fstout "$DIR_OUT/$pop3-$pop2"
# future analyses should rely on 'weighted' Fst (i.e. the ratio of expectations across loci).
realSFS fst stats "$DIR_OUT/$pop3-$pop2.fst.idx" > "$DIR_OUT/$pop3-$pop2.fst.stats"
# summary fst
realSFS fst stats2 "$DIR_OUT/$pop3-$pop2.fst.idx" > "$DIR_OUT/$pop3-$pop2.fstAll"
# windows fst
realSFS fst stats2 "$DIR_OUT/$pop3-$pop2.fst.idx" -win 5000 -step 1000 > "$DIR_OUT/$pop3-$pop2.fstWindows.stats"

# 3-way Fst output (w/ population branch statistic)
realSFS fst index $pop1.saf.idx $pop2.saf.idx $pop3.saf.idx -sfs $pop1-$pop2.sfs -sfs $pop3-$pop1.sfs -sfs $pop3-$pop2.sfs -fstout $pop1-$pop2-$pop3
echo "all done!"
# future analyses should rely on 'weighted' Fst (i.e. the ratio of expectations across loci).
realSFS fst stats "$DIR_OUT/$pop1-$pop2-$pop3.fst.idx" > "$DIR_OUT/$pop1-$pop2-$pop3.fst.stats"
# summary fst
realSFS fst stats2 "$DIR_OUT/$pop1-$pop2-$pop3.fst.idx" > "$DIR_OUT/$pop1-$pop2-$pop3.fstAll"
# windows fst
realSFS fst stats2 "$DIR_OUT/$pop1-$pop2-$pop3.fst.idx" -win 5000 -step 1000 > "$DIR_OUT/$pop1-$pop2-$pop3.fstWindows.stats"


# options
# basic quality filtering for reads
# anc polarizes SFS by reference genome
# underFlowProtect is necessary for large #s of bams
