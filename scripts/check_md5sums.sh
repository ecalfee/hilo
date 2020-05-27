#!/bin/bash
#SBATCH --partition=low2
#SBATCH -D /home/ecalfee/hilo
#SBATCH -J md5_fq
#SBATCH -o /home/ecalfee/hilo/slurm-log/check_md5sums_%A_%a.out
#SBATCH -t 30:00
#SBATCH --mem=4G
#SBATCH --array=1-200

# to run: sbatch check_md5sums.sh
i=$SLURM_ARRAY_TASK_ID
echo "HILO$i checking MD5 sums for fastq files"
md5sum /group/jrigrp6/DanAlignments/HILO"$i"/*.fq.gz > data/HILO_raw_reads/MD5_checksums/March2018/HILO"$i"_MD5.txt
