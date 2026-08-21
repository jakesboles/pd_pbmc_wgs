#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name bowtie2
#SBATCH --nodes 1
#SBATCH --array=1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 32G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1042/Gate_Lab/boles/pd_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

cd /projects/b1042/Gate_Lab/boles/pd_wgs/trimmed_fastqs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "../bowtie_params_id.txt")
R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "../bowtie_params_r1.txt")
R2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "../bowtie_params_r2.txt")

echo ${SAMPLE}

module load bowtie2/2.5.4
module load samtools/1.16.1-gcc-10.4.0

bowtie2 \
-x /projects/p31535/boles/bowtie2_ref_v2/Homo_sapiens_assembly38 \
-p 16 \
-1 ${R1} \
-2 ${R2} \
2> ../bowtie2_reports/${SAMPLE}.bowtie2.log \
| samtools sort -@ 8 -o ../bam/${SAMPLE}.sorted.bam

samtools index ../bam/${SAMPLE}.sorted.bam


