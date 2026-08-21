#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name haplotype_caller
#SBATCH --nodes 1
#SBATCH --array=2-22
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1042/Gate_Lab/boles/pd_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

module load gatk/4.4.0.0

cd /projects/b1042/Gate_Lab/boles/pd_wgs

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "bowtie_params_id.txt")
BAM="bwa_bam/${SAMPLE}.bqsr.bam"

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /projects/p31535/boles/Homo_sapiens_assembly38.fasta \
   -I ${BAM} \
   -O haplotype_caller/${SAMPLE}.output.g.vcf.gz \
   -ERC GVCF