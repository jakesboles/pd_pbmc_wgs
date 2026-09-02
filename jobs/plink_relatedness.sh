#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name plink_relatedness
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 16G
#SBATCH --time 2:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A.log
#SBATCH --verbose

# Estimates pairwise relatedness across the WGS cohort as an additional
# QC check, independent of the WGS<->scATAC identity crosscheck (steps
# 21-24 in CLAUDE.md): this looks for *cryptic relatedness between WGS
# subjects* (duplicate enrollments, unreported family relationships),
# not identity between datasets. Not a SLURM array -- this is a single
# set of pairwise comparisons across the whole 121-sample cohort at once,
# the same shape as gatk_filter_split.sh.
#
# Uses PLINK2's KING-robust kinship estimator (--make-king-table), which
# -- unlike classic IBD/PI_HAT estimation (plink1.9's --genome) -- is
# robust to population stratification, appropriate for a cohort that
# isn't necessarily one homogeneous ancestry group. Kinship coefficients
# are on the KING scale: ~0.5 duplicate/MZ twin, ~0.25 1st-degree,
# ~0.125 2nd-degree, ~0.0625 3rd-degree (halving each step out); the
# conventional midpoint cutoffs for calling a category are ~0.354/0.177/
# 0.0884/0.0442 respectively. No --king-table-filter is set, so the
# output includes every pairwise comparison (~7260 for 121 samples), not
# just flagged/related ones -- useful here since we also want to confirm
# the bulk of pairs cluster near 0, not just spot-check the outliers.
#
# --new-id-max-allele-len 1000 truncate: --set-all-var-ids errors out by
# default on any variant whose REF/ALT allele is longer than its (small)
# built-in cap, which real structural indels in a 121-genome cohort will
# exceed -- 1000 is generous headroom so essentially nothing hits it in
# practice. "truncate" (not the "missing" mode PLINK2 suggests in its own
# error message) is used deliberately: "missing" sets every over-length
# variant's ID to the literal '.', and since --extract later depends on
# these IDs being unique to select specific pruned-in variants, hundreds
# of otherwise-unrelated variants silently sharing one ID would make that
# selection ambiguous. Truncating the allele string in the ID instead
# keeps effectively-unique (if imprecise) IDs for that rare tail.
#
# QC filters below (--maf 0.05 --geno 0.05 --mind 0.1) are standard
# defaults, but run as two SEPARATE, ORDERED plink2 calls -- --geno
# (drop poorly-genotyped variants) strictly before --mind (drop
# poorly-genotyped samples), not combined into one call. This isn't
# cosmetic: PLINK2 computes per-sample missingness against whatever
# variant set is currently loaded, so if --mind runs before --geno has
# removed the worst sites, a relatively small number of bad/low-
# confidence variants (routine in a ~24M-variant joint-genotyped
# multi-sample VCF, where any given sample can easily lack a confident
# call at a site private to other samples) drags every sample's apparent
# missingness rate up -- this is exactly what happened on the first
# attempt here, with --mind computed against the full unfiltered variant
# set and removing all 121 samples. Filtering variants first, then
# computing --mind against that cleaned-up set, is standard practice in
# published GWAS QC protocols for this reason -- these two steps are not
# commutative. Deliberately NOT applying --hwe here: Hardy-Weinberg
# deviation is expected at real sites when a cohort contains related
# individuals -- exactly what this step exists to detect -- so filtering
# on it first would be circular; save --hwe for later analyses that
# assume unrelated samples, after this step has actually established
# that. Restricted to autosomes (--autosome) since chrX-based kinship
# needs per-sample sex, which isn't tracked anywhere in this pipeline,
# and the autosomes alone are more than sufficient for relatedness
# inference.
#
# LD-pruning (--indep-pairwise) before kinship estimation is the
# conventional step in most QC protocols, though KING-robust itself is
# less sensitive to background LD than PCA or PI_HAT are -- treat the
# window/step/r^2 values below as reasonable, adjustable defaults, not
# the only correct choice.

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load plink/2.001

VCF="vqsr/cohort.pass.normalized.vcf.gz"
OUT_DIR="relatedness"

mkdir -p "${OUT_DIR}"

echo "Importing VCF to PLINK2 binary format"

plink2 \
  --threads 8 \
  --vcf "${VCF}" \
  --double-id \
  --max-alleles 2 \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 1000 truncate \
  --autosome \
  --make-pgen \
  --out "${OUT_DIR}/cohort_raw"

echo "Filtering low-quality/rare variants (--maf, --geno)"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_raw" \
  --maf 0.05 \
  --geno 0.05 \
  --make-pgen \
  --out "${OUT_DIR}/cohort_geno"

echo "Filtering poorly-genotyped samples (--mind)"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_geno" \
  --mind 0.1 \
  --make-pgen \
  --out "${OUT_DIR}/cohort_qc"

# --maf/--geno above only ever drop variants, never samples, so the drop
# in sample count between cohort_raw and cohort_qc is attributable
# entirely to --mind. .psam files have one header line (starting with
# '#') followed by one row per sample.
n_before=$(grep -vc '^#' "${OUT_DIR}/cohort_raw.psam")
n_after=$(grep -vc '^#' "${OUT_DIR}/cohort_qc.psam")
echo "--mind 0.1 removed $((n_before - n_after)) of ${n_before} samples (${n_after} remain)"

echo "LD-pruning to an approximately independent marker set"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_qc" \
  --indep-pairwise 200 50 0.1 \
  --out "${OUT_DIR}/cohort_pruned"

echo "Estimating pairwise KING-robust kinship"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_qc" \
  --extract "${OUT_DIR}/cohort_pruned.prune.in" \
  --make-king-table \
  --out "${OUT_DIR}/cohort_king"

echo "Done"
