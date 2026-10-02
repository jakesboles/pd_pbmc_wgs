#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name crosscheck_fingerprints
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 32G
#SBATCH --time 8:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%j.log
#SBATCH --verbose

# Verifies donor identity between the filtered WGS cohort VCF and every
# matched scATAC-seq (Cell Ranger ARC) possorted BAM. Reworked from a
# per-pair SLURM array (one task per presumed WGS<->ATAC match, 121 tasks)
# into a single all-pairs comparison, specifically to catch sample-label
# swaps that the old per-pair design could not.
#
# WHY THE OLD DESIGN COULD NOT DETECT A SWAP: each task subsetted --INPUT
# down to exactly one VCF sample (SelectVariants -sn <wgs_sample>) and
# --SECOND_INPUT down to exactly that one sample's presumed-matching BAM,
# then asked "do these two match." If, say, JSB100-1's true genetic match
# was actually the BAM labeled for JSB100-2 (a real label swap somewhere
# upstream in sample handling/demographics), that task would just report
# EXPECTED_MISMATCH -- it never had a chance to notice the real match
# sitting one row away, because that BAM was never part of the
# comparison at all. 121 independent one-pair tasks can only confirm or
# deny each presumed pairing in isolation; they structurally cannot
# surface a swap between two samples.
#
# This version passes the FULL cohort VCF (all 121 samples, one --INPUT)
# against ALL 121 reordered/subset ATAC BAMs (repeated --SECOND_INPUT
# arguments, built below from params/crosscheck_sample_map.txt) in one
# job, with CROSSCHECK_MODE CHECK_ALL_OTHERS instead of the default
# CHECK_SAME_SAMPLE. CHECK_ALL_OTHERS does everything CHECK_SAME_SAMPLE
# does (confirms each sample's presumed match) AND additionally confirms
# that every sample does NOT unexpectedly match any OTHER sample --
# i.e. it computes and reports the full 121x121 LOD matrix, not just the
# 121 presumed-matched pairs. Each row of the output is one
# VCF-sample x BAM-sample comparison; RESULT is EXPECTED_MATCH/
# EXPECTED_MISMATCH when the presumed pairing (per INPUT_SAMPLE_MAP) came
# back as expected, or UNEXPECTED_MATCH/UNEXPECTED_MISMATCH when it
# didn't -- an UNEXPECTED_MATCH between two different presumed samples is
# exactly the sample-swap signature to go looking for.
#
# --INPUT_SAMPLE_MAP is still needed, and is still the full 121-row map
# as-is: with both sides now carrying all 121 samples, every VCF sample
# has a real counterpart to compare against on the ATAC side (matched by
# scATAC SM tag, not the "JSB"-prefixed WGS name), so the old per-task
# "sample X is missing from RIGHT group" log noise doesn't recur here --
# that noise was a symptom of a 121-vs-1 size mismatch between INPUT and
# SECOND_INPUT, which no longer exists once both sides carry the full
# cohort. No per-sample VCF subsetting (SelectVariants -sn) is needed any
# more either, for the same reason: the point now is to compare every
# sample against every other sample, not isolate one pairing at a time.
#
# --SECOND_INPUT is built as an explicit, repeated argument per BAM
# (looped below) rather than relying on this GATK build's .list-file
# support for BAM inputs, which isn't confirmed on this cluster -- an
# explicitly generated argument list is unambiguous either way.
#
# Heavier than the old per-pair tasks (121x121 ~= 14,641 pairwise LOD
# computations in one job instead of 121 jobs of 1 each), so bumped to
# 32G/8h from the old 16G/2h. Not a SLURM array any more -- there's only
# one job to submit. Still depends on subset_reorder_atac_bams.sh having
# already produced crosscheck/atac_subset/<wgs_sample>.subset.reordered.bam
# for every sample in params/crosscheck_sample_map.txt.

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load gatk/4.4.0.0

HAPLOTYPE_MAP="/projects/p31535/boles/Homo_sapiens_assembly38.haplotype_database.txt"
VCF="vqsr/cohort.pass.normalized.vcf.gz"
SAMPLE_MAP="params/crosscheck_sample_map.txt"

mkdir -p crosscheck

SECOND_INPUT_ARGS=()
n_bams=0
while IFS=$'\t' read -r wgs_sample atac_sample; do
  reordered_bam="crosscheck/atac_subset/${wgs_sample}.subset.reordered.bam"
  if [[ ! -f "${reordered_bam}" ]]; then
    echo "ERROR: missing reordered BAM for ${wgs_sample}: ${reordered_bam}" >&2
    echo "Run subset_reorder_atac_bams.sh for this sample first." >&2
    exit 1
  fi
  SECOND_INPUT_ARGS+=(--SECOND_INPUT "${reordered_bam}")
  n_bams=$((n_bams + 1))
done < "${SAMPLE_MAP}"

echo "Comparing all samples in ${VCF} against ${n_bams} reordered ATAC BAMs (full all-pairs matrix)"

# CROSSCHECK_BY SAMPLE (default is READGROUP) so each sample is compared
# once as a whole, not once per read group. EXIT_CODE_WHEN_MISMATCH is 0
# because a genotype mismatch -- including an UNEXPECTED_MATCH flagging a
# real sample swap -- is an expected possible QC finding to review in the
# output metrics, not a pipeline failure; EXIT_CODE_WHEN_NO_VALID_CHECKS
# is left at its default so a real misconfiguration (e.g. no overlapping
# fingerprinting sites) still fails the job loudly.
gatk CrosscheckFingerprints \
  --INPUT "${VCF}" \
  "${SECOND_INPUT_ARGS[@]}" \
  --INPUT_SAMPLE_MAP "${SAMPLE_MAP}" \
  --HAPLOTYPE_MAP "${HAPLOTYPE_MAP}" \
  --CROSSCHECK_MODE CHECK_ALL_OTHERS \
  --CROSSCHECK_BY SAMPLE \
  --EXIT_CODE_WHEN_MISMATCH 0 \
  --OUTPUT "crosscheck/cohort_all_pairs.crosscheck_metrics"

echo "Done -- review crosscheck/cohort_all_pairs.crosscheck_metrics for any UNEXPECTED_MATCH/UNEXPECTED_MISMATCH RESULT"
