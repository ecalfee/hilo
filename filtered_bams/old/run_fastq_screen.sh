#!/bin/bash
#SBATCH --partition=med2
#SBATCH -D /home/ecalfee/hilo/filtered_bams
#SBATCH -J fastq_screen
#SBATCH -o /home/ecalfee/hilo/slurm-log/fastq_screen_%A_%a.out
#SBATCH -t 4:00:00
#SBATCH --mem=16G
#SBATCH --array=0
#SBATCH --export=DIR_IN=../data/HILO_raw_reads,PREFIX_LIST=Jan2019

# START WITH ZERO FOR ARRAY INDEXES (!)

# this script takes in fastq files with raw paired reads
# and aligns reads to APGv4 reference with all chromosome and contigs
# sample ID, fastq file 1, fastq file 2, output bam directory
# to run: sbatch --array=1 --export=DIR_IN=../data/HILO_raw_reads/TEST,PREFIX_LIST=Jan2019 map_and_filter_reads.sh

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# load modules
module load bio3
module load fastq_screen

REF="../data/refMaize/Zea_mays.B73_RefGen_v4.dna.toplevel.fa"
i=$SLURM_ARRAY_TASK_ID
ids=($(cat ${DIR_IN}/${PREFIX_LIST}_IDs.list))
ID=${ids["$i"]}
libraries=($(cat ${DIR_IN}/${PREFIX_LIST}_libraries.list))
LIBRARY=${libraries["$i"]}
lanes=($(cat ${DIR_IN}/${PREFIX_LIST}_lanes.list))
LANE=${lanes["$i"]}
FASTQ1=${DIR_IN}/${LANE}/${ID}_1.fq.gz
FASTQ2=${DIR_IN}/${LANE}/${ID}_2.fq.gz
DIR_OUT=results/${LANE} # results directory for final and intermediate bam files

echo "all done!"
fastq_screen --outdir=fastQC --aligner BOWTIE2 --conf /home/ecalfee/hilo/filtered_bams/fastq_screen/fastq_screen_maize.conf --outdir fastq_screen/results ../../../../group/jrigrp10/ErinHilo/data/MTPJRI01/HILO300_S90_L002_R1_001.fastq.gz
# options
# default is mapping a subset of 100,000 reads
# configuration file specifies genomes to map to:
# human, (human) mtdna, mouse, rat, ecoli, lambda, drosophila, phix, adapters, vectors, rRNA, worm, yeast, arabidopsis, maize
