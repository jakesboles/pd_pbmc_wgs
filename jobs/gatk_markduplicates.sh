#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name gatk_markduplicates
#SBATCH --nodes 1
#SBATCH --array=59,99
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1042/Gate_Lab/boles/pd_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

cd /projects/b1042/Gate_Lab/boles/pd_wgs

module load samtools/1.16.1-gcc-10.4.0
module load gatk/4.4.0.0

PARAMS_FILE="bowtie_params_id.txt"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f1 -d,)

file=$(ls "$PWD"/bwa_bam/* | grep "${sample}[.]" | grep ".merged.sorted")

echo ${sample}
echo ${file}

gatk MarkDuplicates \
-I ${file} \
-O bwa_bam/${sample}.markdup.bam \
-M gatk_reports/${sample}.duplicate_metrics.txt

echo "Marked duplicates successfully"

samtools index -@ 16 bwa_bam/${sample}.markdup.bam

echo "Indexed duplicate .bam successfully"