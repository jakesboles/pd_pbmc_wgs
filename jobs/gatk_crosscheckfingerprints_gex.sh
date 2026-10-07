#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name crosscheck_fingerprints_gex
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 32G
#SBATCH --time 8:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%j.log
#SBATCH --verbose

# Runs under b1169/b1169 (not b1042/genomics), matching
# gatk_crosscheckfingerprints.sh -- the ATAC version of this all-pairs job
# needed the extra memory/tmp headroom in practice, and this job is the
# same shape (one job, full cohort VCF against all matched BAMs at once),
# so it's given the same allocation up front rather than waiting to hit
# the same wall.
#
# Verifies donor identity between the filtered WGS cohort VCF and every
# matched GEX (scRNA) BAM -- an independent additional sanity check on
# WGS<->multiome identity, alongside the ATAC-based check in
# gatk_crosscheckfingerprints.sh. Same all-pairs design as that script
# (see its header comment for the full rationale on why an all-pairs
# comparison, not a per-presumed-pair one, is what's needed to actually
# rule out a sample label swap): one job, --INPUT is the full cohort VCF
# (all 121 samples), --SECOND_INPUT is repeated once per sample for all
# matched reordered/subset GEX BAMs
# (crosscheck/gex_subset/<wgs_sample>.subset.reordered.bam, from
# subset_reorder_gex_bams.sh), --INPUT_SAMPLE_MAP renames each WGS VCF
# sample to its GEX RG SM tag, CROSSCHECK_MODE CHECK_ALL_OTHERS computes
# the full LOD/RESULT matrix so any UNEXPECTED_MATCH (the sample-swap
# signature) would be visible.
#
# CAVEAT SPECIFIC TO THE GEX COMPARISON: the haplotype map's SNP sites
# were chosen for roughly uniform genome-wide coverage, which is a good
# assumption for WGS and for ATAC fragments (open-chromatin-biased but
# still genome-wide), but NOT for RNA-seq -- GEX reads only cover
# transcribed, predominantly exonic sequence, and read depth at any
# given site tracks that gene's expression level, not genomic position.
# So expect meaningfully FEWER informative (covered) haplotype-map sites
# per sample here than in the ATAC comparison, and correspondingly
# weaker (though still decisive for a true match, assuming enough
# covered sites survive) LOD scores -- this is an expected property of
# comparing against RNA-seq data, not a bug. Genes with very low/no
# expression in PBMCs contribute nothing here, independent of any
# sample-identity question.
#
# No special RNA-seq preprocessing (e.g. SplitNCigarReads) is applied --
# CrosscheckFingerprints does its own pileup-based genotype-likelihood
# comparison at the given sites rather than calling variants, so spliced
# reads are not expected to need special handling here, consistent with
# how the ATAC BAMs were used directly (aside from the subset/reorder
# step both share).

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load gatk/4.4.0.0

HAPLOTYPE_MAP="/projects/p31535/boles/Homo_sapiens_assembly38.haplotype_database.txt"
VCF="vqsr/cohort.pass.normalized.vcf.gz"
SAMPLE_MAP="params/crosscheck_gex_sample_map.txt"

mkdir -p crosscheck

SECOND_INPUT_ARGS=()
n_bams=0
while IFS=$'\t' read -r wgs_sample gex_sample; do
  reordered_bam="crosscheck/gex_subset/${wgs_sample}.subset.reordered.bam"
  if [[ ! -f "${reordered_bam}" ]]; then
    echo "ERROR: missing reordered BAM for ${wgs_sample}: ${reordered_bam}" >&2
    echo "Run subset_reorder_gex_bams.sh for this sample first." >&2
    exit 1
  fi
  SECOND_INPUT_ARGS+=(--SECOND_INPUT "${reordered_bam}")
  n_bams=$((n_bams + 1))
done < "${SAMPLE_MAP}"

echo "Comparing all samples in ${VCF} against ${n_bams} reordered GEX BAMs (full all-pairs matrix)"

gatk CrosscheckFingerprints \
  --INPUT "${VCF}" \
  "${SECOND_INPUT_ARGS[@]}" \
  --INPUT_SAMPLE_MAP "${SAMPLE_MAP}" \
  --HAPLOTYPE_MAP "${HAPLOTYPE_MAP}" \
  --CROSSCHECK_MODE CHECK_ALL_OTHERS \
  --CROSSCHECK_BY SAMPLE \
  --EXIT_CODE_WHEN_MISMATCH 0 \
  --OUTPUT "crosscheck/cohort_all_pairs_gex.crosscheck_metrics"

echo "Done -- review crosscheck/cohort_all_pairs_gex.crosscheck_metrics for any UNEXPECTED_MATCH/UNEXPECTED_MISMATCH RESULT"
