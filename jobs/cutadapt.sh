#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name cutadapt
#SBATCH --nodes 1
#SBATCH --array=1-902
#SBATCH --ntasks-per-node 16
#SBATCH --mem 32G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

module load cutadapt/4.2

cd /projects/b1169/boles/pd_pbmc_wgs

PARAMS_FILE="cutadapt_params.txt"

R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f1 -d,)
R2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f2 -d,)
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f3 -d,)

cutadapt \
  -j 16 \
  -a CTGTCTCTTATACACATCT \
  -A CTGTCTCTTATACACATCT \
  --nextseq-trim 20 \
  --minimum-length 20 \
  -o trimmed_fastqs/${R1} \
  -p trimmed_fastqs/${R2} \
  "fastq/${R1}" "fastq/${R2}" \
  > cutadapt_logs/${sample}.cutadapt.log
