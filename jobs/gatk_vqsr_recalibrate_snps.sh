#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name vqsr_snp_recal
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 128G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1042/Gate_Lab/boles/pd_wgs/logs/%x_%A.log
#SBATCH --verbose

cd /projects/b1042/Gate_Lab/boles/pd_wgs

module load samtools/1.16.1-gcc-10.4.0
module load gatk/4.4.0.0
module load R/4.4.0

# Hardcode Java 17 explicitly — don't trust PATH resolution
JAVA17_BIN=/software/java/jdk-17.0.2+8/bin/java    # <-- fill in from module show output

${JAVA17_BIN} -version   # confirm this reports 17.0.2 before proceeding

# ---- Step 1: index the input VCF ----
echo "Indexing input VCF"
${JAVA17_BIN} -jar ${GATK_DIR}/gatk-package-4.4.0.0-local.jar IndexFeatureFile \
   -I gathered_genotyped_gvcf/genotyped_cohort.raw.vcf.gz

# ---- Step 2: VariantRecalibratory (SNP model) ----
echo "Running VariantRecalibrator with SNPs"
${JAVA17_BIN} -Xmx32g -jar ${GATK_DIR}/gatk-package-4.4.0.0-local.jar VariantRecalibrator \
   -R /projects/p31535/boles/Homo_sapiens_assembly38.fasta \
   -V gathered_genotyped_gvcf/genotyped_cohort.raw.vcf.gz \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /projects/p31535/boles/hapmap_3.3.hg38.vcf.gz \
   --resource:omni,known=false,training=true,truth=true,prior=12.0 /projects/p31535/boles/1000G_omni2.5.hg38.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 /projects/p31535/boles/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 /projects/p31535/boles/dbsnp_146.hg38.vcf.gz \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
   -mode SNP \
   --tranches-file vqsr/cohort.snps.tranches \
   --rscript-file vqsr/cohort.snps.plots.R \
   -O vqsr/cohort.snps.recal