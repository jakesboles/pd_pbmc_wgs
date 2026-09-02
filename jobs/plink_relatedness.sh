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
# QC filters below (--maf 0.05 --geno 0.05 --mind 0.1) are standard
# defaults. Deliberately NOT applying --hwe here: Hardy-Weinberg
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
  --autosome \
  --make-pgen \
  --out "${OUT_DIR}/cohort_raw"

echo "Applying variant/sample QC filters"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_raw" \
  --maf 0.05 \
  --geno 0.05 \
  --mind 0.1 \
  --make-pgen \
  --out "${OUT_DIR}/cohort_qc"

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
