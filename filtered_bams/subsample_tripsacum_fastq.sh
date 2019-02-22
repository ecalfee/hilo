#!/bin/bash
#SBATCH --partition=bigmemm
#SBATCH -D /home/ecalfee/hilo/data/tripsacum
#SBATCH -J sampleTrip
#SBATCH -o /home/ecalfee/hilo/slurm-log/subsampleTripsacumFastQ_%A_%a.out
#SBATCH -t 2-00:00:00
#SBATCH --mem=8G
#SBATCH --array=1-2
#SBATCH --export=FRACTION=0.5,SRR=SRR7758238

# this script takes a subsample of fraction $FRACTION reads from tripsacum

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio # for seqtk

echo "working directory: ${PWD}" # print current working directory

i=$SLURM_ARRAY_TASK_ID
SEED=7891

# make directory for output
mkdir -p "$SRR"

# sampling reads
echo "sampling reads, fraction = $FRACTION"
seqtk sample -s$SEED "$SRR"_"$i".fastq.gz $FRACTION | gzip > "$SRR"/TRIP_"$i".fq.gz

echo "all done!"
# note: random seed must be consistent between _1 and _2 paired read files

