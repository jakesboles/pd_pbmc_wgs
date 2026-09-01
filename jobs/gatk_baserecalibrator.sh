#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name gatk_recalibrate
#SBATCH --nodes 1
#SBATCH --array=1-121
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

cd /projects/b1169/boles/pd_pbmc_wgs

module load samtools/1.16.1-gcc-10.4.0
module load gatk/4.4.0.0

PARAMS_FILE="bowtie_params_id.txt"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f1 -d,)

echo ${sample}

echo "Recalibrating reads"

gatk BaseRecalibrator \
    -R /projects/p31535/boles/Homo_sapiens_assembly38.fasta \
    -I bwa_bam/${sample}.markdup.bam \
    --known-sites /projects/p31535/boles/dbsnp_146.hg38.vcf.gz \
    --known-sites /projects/p31535/boles/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O gatk_reports/${sample}.recal.table
    
echo "Applying recalibration"
    
gatk ApplyBQSR \
    -R /projects/p31535/boles/Homo_sapiens_assembly38.fasta \
    -I bwa_bam/${sample}.markdup.bam \
    --bqsr-recal-file gatk_reports/${sample}.recal.table \
    -O bwa_bam/${sample}.bqsr.bam
