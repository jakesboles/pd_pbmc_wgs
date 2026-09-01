#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name gather_genotype_gvcfs
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 128G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1042/Gate_Lab/boles/pd_wgs/logs/%x_%A.log
#SBATCH --verbose

module load gatk/4.4.0.0

cd /projects/b1042/Gate_Lab/boles/pd_wgs

gatk --java-options "-Djava.io.tmpdir=/projects/b1169/boles/pd_wgs/temp" GatherVcfs \
  -I genotyped_gvcfs/chr1.vcf.gz \
  -I genotyped_gvcfs/chr2.vcf.gz \
  -I genotyped_gvcfs/chr3.vcf.gz \
  -I genotyped_gvcfs/chr4.vcf.gz \
  -I genotyped_gvcfs/chr5.vcf.gz \
  -I genotyped_gvcfs/chr6.vcf.gz \
  -I genotyped_gvcfs/chr7.vcf.gz \
  -I genotyped_gvcfs/chr8.vcf.gz \
  -I genotyped_gvcfs/chr9.vcf.gz \
  -I genotyped_gvcfs/chr10.vcf.gz \
  -I genotyped_gvcfs/chr11.vcf.gz \
  -I genotyped_gvcfs/chr12.vcf.gz \
  -I genotyped_gvcfs/chr13.vcf.gz \
  -I genotyped_gvcfs/chr14.vcf.gz \
  -I genotyped_gvcfs/chr15.vcf.gz \
  -I genotyped_gvcfs/chr16.vcf.gz \
  -I genotyped_gvcfs/chr17.vcf.gz \
  -I genotyped_gvcfs/chr18.vcf.gz \
  -I genotyped_gvcfs/chr19.vcf.gz \
  -I genotyped_gvcfs/chr20.vcf.gz \
  -I genotyped_gvcfs/chr21.vcf.gz \
  -I genotyped_gvcfs/chr22.vcf.gz \
  -I genotyped_gvcfs/chrX.vcf.gz \
  -O gathered_genotyped_gvcf/genotyped_cohort.raw.vcf.gz \
  --TMP_DIR /projects/b1169/boles/pd_wgs/temp \
  --MAX_RECORDS_IN_RAM 500000000
