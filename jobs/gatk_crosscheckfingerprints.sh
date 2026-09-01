#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name crosscheck_fingerprints
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A.log
#SBATCH --verbose

# Verifies donor identity between the filtered WGS cohort VCF and each
# matched scATAC-seq (Cell Ranger ARC) possorted BAM, using GATK's
# haplotype-map-based fingerprinting. Requires crosscheck_sample_map.txt
# and crosscheck_atac_bams.txt from make_crosscheck_params.sh to exist.

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load gatk/4.4.0.0

HAPLOTYPE_MAP="/projects/p31535/boles/Homo_sapiens_assembly38.haplotype_database.txt"
VCF="vqsr/cohort.pass.normalized.vcf.gz"

mkdir -p crosscheck

# One --SECOND_INPUT flag per ATAC BAM in the crosswalk. Passed this way
# (rather than pointing GATK at the .txt file directly) so behavior
# doesn't depend on whether this GATK version still honors Picard's
# ".list"-extension file-of-files convention.
SECOND_INPUT_ARGS=()
while read -r bam; do
  SECOND_INPUT_ARGS+=(--SECOND_INPUT "${bam}")
done < crosscheck_atac_bams.txt

# INPUT_SAMPLE_MAP renames each WGS VCF sample (e.g. JSB100-1) to the
# scATAC sample name recorded in the matching BAM's RG SM tag (e.g.
# 100-1), so CHECK_SAME_SAMPLE can pair them up despite the "JSB" prefix
# mismatch. CROSSCHECK_BY is forced to SAMPLE (default is READGROUP) so
# each WGS sample is compared once against its whole matched BAM rather
# than once per read group. EXIT_CODE_WHEN_MISMATCH is set to 0 because a
# genotype mismatch here is an expected possible QC finding (e.g. a
# sample swap) to review in the output metrics, not a pipeline failure —
# EXIT_CODE_WHEN_NO_VALID_CHECKS is left at its default so a real
# misconfiguration (e.g. no overlapping fingerprinting sites) still fails
# the job loudly.
gatk CrosscheckFingerprints \
  --INPUT "${VCF}" \
  "${SECOND_INPUT_ARGS[@]}" \
  --INPUT_SAMPLE_MAP crosscheck_sample_map.txt \
  --HAPLOTYPE_MAP "${HAPLOTYPE_MAP}" \
  --CROSSCHECK_MODE CHECK_SAME_SAMPLE \
  --CROSSCHECK_BY SAMPLE \
  --NUM_THREADS 16 \
  --EXIT_CODE_WHEN_MISMATCH 0 \
  --OUTPUT crosscheck/cohort_vs_atac.crosscheck_metrics
