#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name vqsr_indel_apply
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 128G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A.log
#SBATCH --verbose

cd /projects/b1169/boles/pd_pbmc_wgs

module load samtools/1.16.1-gcc-10.4.0
module load gatk/4.4.0.0
module load R/4.4.0

# Hardcode Java 17 explicitly — don't trust PATH resolution
JAVA17_BIN=/software/java/jdk-17.0.2+8/bin/java    # <-- fill in from module show output

${JAVA17_BIN} -version   # confirm this reports 17.0.2 before proceeding

OUT_DIR="vqsr"
RESOURCE_DIR="/projects/p31535/boles"

# ---- Step 5: ApplyVQSR (Indel) ----
echo "Applying indel VQSR model"
${JAVA17_BIN} -Xmx12g -jar ${GATK_DIR}/gatk-package-4.4.0.0-local.jar ApplyVQSR \
   -R ${RESOURCE_DIR}/Homo_sapiens_assembly38.fasta \
   -V ${OUT_DIR}/cohort.snps_recalibrated.vcf.gz \
   --recal-file ${OUT_DIR}/cohort.indels.recal \
   --tranches-file ${OUT_DIR}/cohort.indels.tranches \
   --truth-sensitivity-filter-level 99.0 \
   -mode INDEL \
   -O ${OUT_DIR}/cohort.recalibrated.vcf.gz