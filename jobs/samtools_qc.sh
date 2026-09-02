#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name samtools_qc
#SBATCH --nodes 1
#SBATCH --array=59,99
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A.log
#SBATCH --verbose

module load samtools/1.16.1-gcc-10.4.0

cd /projects/b1169/boles/pd_pbmc_wgs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "params/bowtie_params_id.txt")

echo ${SAMPLE}

echo "Running flagstat"

samtools flagstat -@ 16 bwa_bam/${SAMPLE}.merged.sorted.bam \
    > samtools_reports/${SAMPLE}.flagstat.txt

echo "Running stats"

samtools stats -@ 16 bwa_bam/${SAMPLE}.merged.sorted.bam \
    > samtools_reports/${SAMPLE}.stats.txt

echo "Running idxstats"

samtools idxstats bwa_bam/${SAMPLE}.merged.sorted.bam \
    > samtools_reports/${SAMPLE}.idxstats.txt

echo "Running coverage"

samtools coverage bwa_bam/${SAMPLE}.merged.sorted.bam \
    > samtools_reports/${SAMPLE}.coverage.txt
    
echo "Making duplicates report"
    
samtools flagstat -@ 16 bwa_bam/${SAMPLE}.markdup.bam \
    > samtools_reports/${SAMPLE}.markdup.flagstat.txt
    
samtools stats -@ 16 bwa_bam/${SAMPLE}.markdup.bam \
    > samtools_reports/${SAMPLE}.markdup.stats.txt
samtools idxstats bwa_bam/${SAMPLE}.markdup.bam \
    > samtools_reports/${SAMPLE}.markdup.idxstats.txt
samtools coverage bwa_bam/${SAMPLE}.markdup.bam \
    > samtools_reports/${SAMPLE}.markdup.coverage.txt
    
# echo "Getting depth"
#     
# samtools depth bam/${SAMPLE}.sorted.bam \
# | awk '{sum+=$3; n++} END {print sum/n}'
