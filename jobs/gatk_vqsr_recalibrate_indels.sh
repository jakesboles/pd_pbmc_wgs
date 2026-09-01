#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name vqsr_indel_recal
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

# ---- Step 4: VariantRecalibrator (Indel model) ----
echo "Running VariantRecalibrator with indels"
${JAVA17_BIN} -Xmx24g -jar ${GATK_DIR}/gatk-package-4.4.0.0-local.jar VariantRecalibrator \
   -R ${RESOURCE_DIR}/Homo_sapiens_assembly38.fasta \
   -V ${OUT_DIR}/cohort.snps_recalibrated.vcf.gz \
   --resource:mills,known=false,training=true,truth=true,prior=12.0 ${RESOURCE_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${RESOURCE_DIR}/dbsnp_146.hg38.vcf.gz \
   -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an DP \
   -mode INDEL \
   --tranches-file ${OUT_DIR}/cohort.indels.tranches \
   --rscript-file ${OUT_DIR}/cohort.indels.plots.R \
   -O ${OUT_DIR}/cohort.indels.recal