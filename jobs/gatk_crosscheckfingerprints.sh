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
# and crosscheck_atac_bams.txt from make_crosscheck_params.sh, AND the
# reordered subset BAMs from subset_reorder_atac_bams.sh, to already
# exist -- with matching row order/count (one array task per line/pair).
#
# One task per sample pair, rather than one job comparing the VCF against
# all matched BAMs at once, so a problem with any single comparison (e.g.
# a real genotype mismatch/sample swap, or a task-level error) only fails
# that one task instead of the whole cohort.
#
# SECOND_INPUT points at the reordered subset BAM
# (crosscheck/atac_subset/<wgs_sample>.subset.reordered.bam), not the raw
# atac_possorted_bam.bam -- the raw ATAC BAMs list contigs in a different
# order than the WGS VCF/haplotype map (same contigs, alphabetical vs.
# numeric), which fails CrosscheckFingerprints's strict sequence-
# dictionary check for every sample. See subset_reorder_atac_bams.sh for
# the fix (subset to fingerprinting sites, then reorder to match).
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
reordered_bam="crosscheck/atac_subset/${wgs_sample}.subset.reordered.bam"

echo "${wgs_sample}"
echo "${bam}"
echo "${reordered_bam}"

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
  --SECOND_INPUT "${reordered_bam}" \
  --INPUT_SAMPLE_MAP crosscheck_sample_map.txt \
  --HAPLOTYPE_MAP "${HAPLOTYPE_MAP}" \
  --CROSSCHECK_MODE CHECK_SAME_SAMPLE \
  --CROSSCHECK_BY SAMPLE \
  --EXIT_CODE_WHEN_MISMATCH 0 \
  --OUTPUT "crosscheck/${wgs_sample}.crosscheck_metrics"
