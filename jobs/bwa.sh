#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name bwa
#SBATCH --nodes 1
#SBATCH --array=1-902
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

cd /projects/b1169/boles/pd_pbmc_wgs

module load samtools/1.16.1-gcc-10.4.0

PARAMS_FILE="params/bwa_params.txt"

R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f1 -d,)
R2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f2 -d,)
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f3 -d,)
lane=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f4 -d,)
id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f5 -d,)

library="${id}_lib1"
rgid="${id}_${lane}"

header="@RG\tID:${rgid}\tSM:${id}\tLB:${library}\tPL:ILLUMINA\tPU:${lane}"

# echo "${header}"

/projects/p31535/boles/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem \
  -t 16 \
  -R ${header} \
  /projects/p31535/boles/Homo_sapiens_assembly38.fasta \
  "trimmed_fastqs/${R1}" \
  "trimmed_fastqs/${R2}" \
  2> bwa_reports/${sample}.bwa.log \
| samtools sort -@ 16 -o "bwa_bam/${sample}.sorted.bam"

samtools index bwa_bam/${sample}.sorted.bam
