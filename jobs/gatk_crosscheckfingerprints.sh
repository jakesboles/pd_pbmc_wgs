#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name crosscheck_fingerprints
#SBATCH --nodes 1
#SBATCH --array=1-121
#SBATCH --ntasks-per-node 4
#SBATCH --mem 16G
#SBATCH --time 2:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

# Verifies donor identity between the filtered WGS cohort VCF and each
# matched scATAC-seq (Cell Ranger ARC) possorted BAM, using GATK's
# haplotype-map-based fingerprinting. Requires crosscheck_sample_map.txt
# and crosscheck_atac_bams.txt from make_crosscheck_params.sh to exist,
# with matching row order/count (one array task per line/pair).
#
# One task per sample pair, rather than one job comparing the VCF against
# all matched BAMs at once: CrosscheckFingerprints requires every file it
# fingerprints together to share an identical sequence dictionary, and at
# least one ATAC BAM in this cohort was evidently processed against a
# Cell Ranger ARC reference build with a different contig order than the
# rest (a real cross-batch reference inconsistency, not a bug here) --
# batching all 121 BAMs into one CrosscheckFingerprints call meant that
# one mismatched BAM's dictionary error killed the whole job instead of
# just that one comparison. Splitting into an array isolates the failure
# to whichever sample(s) it actually affects; any task that errors out
# (as opposed to completing and reporting an EXPECTED_MISMATCH row) means
# that sample's ATAC BAM needs separate follow-up -- check its Cell
# Ranger ARC reference package/version against the others before trusting
# its multiome data downstream.
#
# NOTE: --array bounds above must match `wc -l crosscheck_atac_bams.txt`
# -- update both if the crosswalk is regenerated with a different sample
# count.

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load gatk/4.4.0.0

HAPLOTYPE_MAP="/projects/p31535/boles/Homo_sapiens_assembly38.haplotype_database.txt"
VCF="vqsr/cohort.pass.normalized.vcf.gz"

mkdir -p crosscheck

wgs_sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" crosscheck_sample_map.txt | cut -f1)
bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" crosscheck_atac_bams.txt)

echo "${wgs_sample}"
echo "${bam}"

# INPUT_SAMPLE_MAP renames each WGS VCF sample (e.g. JSB100-1) to the
# scATAC sample name recorded in the matching BAM's RG SM tag (e.g.
# 100-1), so CHECK_SAME_SAMPLE can pair them up despite the "JSB" prefix
# mismatch -- passing the full map every task is harmless, since only the
# one WGS/ATAC pair present in both this task's INPUT and SECOND_INPUT
# actually gets compared. CROSSCHECK_BY is forced to SAMPLE (default is
# READGROUP) so the WGS sample is compared once against the whole matched
# BAM rather than once per read group. EXIT_CODE_WHEN_MISMATCH is set to
# 0 because a genotype mismatch here is an expected possible QC finding
# (e.g. a sample swap) to review in the output metrics, not a pipeline
# failure -- EXIT_CODE_WHEN_NO_VALID_CHECKS is left at its default so a
# real misconfiguration (e.g. no overlapping fingerprinting sites) still
# fails the task loudly.
gatk CrosscheckFingerprints \
  --INPUT "${VCF}" \
  --SECOND_INPUT "${bam}" \
  --INPUT_SAMPLE_MAP crosscheck_sample_map.txt \
  --HAPLOTYPE_MAP "${HAPLOTYPE_MAP}" \
  --CROSSCHECK_MODE CHECK_SAME_SAMPLE \
  --CROSSCHECK_BY SAMPLE \
  --EXIT_CODE_WHEN_MISMATCH 0 \
  --OUTPUT "crosscheck/${wgs_sample}.crosscheck_metrics"
