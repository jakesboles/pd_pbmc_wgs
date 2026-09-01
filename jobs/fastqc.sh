#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name fastqc
#SBATCH --nodes 1
#SBATCH --array=1-1804
#SBATCH --ntasks-per-node 16
#SBATCH --mem 32G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A.log
#SBATCH --verbose

module load fastqc/0.12.0

cd /projects/b1169/boles/pd_pbmc_wgs

FILE=$(ls fastq | sed -n "${SLURM_ARRAY_TASK_ID}p")

fastqc fastq/${FILE} \
-o fastqc_reports \
--memory 3200 \
--threads 16 \
--noextract
