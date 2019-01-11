#!/bin/bash -l
#SBATCH --partition=med
#SBATCH -D /home/ecalfee/hilo/data
#SBATCH -J out2DSFS
#SBATCH -o /home/ecalfee/hilo/slurm-log/outliers2DSFS_%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --mem=1G
#SBATCH -n 6

# note: must set the directory to print output and find regions.txt file
# and must set array as well when calling the script:
# to run:
# (where 8 is the number of lines in hap_groups.list)
# sbatch --array=1-28 --export="DIR=outliers/chr4/inv4m,ALL" calcFST_2DSFSAngsd_outliers.sh

i=$(($SLURM_ARRAY_TASK_ID-1)) # because 0 indexed
LIST_OF_PAIR1=($(cat $DIR/hap_groups_pair1.list))
POP1=${LIST_OF_PAIR1[{$i}]}
LIST_OF_PAIR2=($(cat $DIR/hap_groups_pair2.files))
POP2=${LIST_OF_PAIR2[{$i}]}

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

mkdir -p "$DIR/2DSFS"
mkdir -p "$DIR/FST"

# calculate pairwise 2D SFS
echo "calculating pairwise 2D SFS for pop $POP1 and pop $POP2 region $DIR"
realSFS "$DIR/SAF/$POP1.saf.idx" "$DIR/SAF/$POP2.saf.idx" -P 6 > "$DIR/SFS/$POP1-$POP2.sfs"

echo "now calculating FST"
# calculate FST (-whichFst 1 indicates Hudson/Bhatia2013 estimator for Fst)
realSFS fst index "$DIR/SAF/$POP1.saf.idx" "$DIR/SAF/$POP2.saf.idx" \
-sfs "$DIR/SFS/$POP1-$POP2.sfs" \
-whichFst 1 \
-fstout "$DIR/FST/$POP1-$POP2"
# future analyses should rely on 'weighted' Fst (i.e. the ratio of expectations across loci).
realSFS fst stats "$DIR/FST/$POP1-$POP2.fst.idx" > "$DIR/FST/$POP1-$POP2.fst.stats"

echo "all done!"
