#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name plink_sex_check_prep
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 16G
#SBATCH --time 2:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A.log
#SBATCH --verbose

# Prep step for the X-heterozygosity biological-sex check. Prep-only,
# same split as jobs/genesis_pcrelate_prep.sh/r_scripts/genesis_pcrelate.R
# -- this job does the PLINK-side data prep; r_scripts/sex_check_viz.R
# does the actual per-sample sex-call computation.
#
# CONFIRMED against this cluster's actual `plink2 --help` output (full
# top-level listing, not the generic web docs): this build -- self-
# reports as PLINK v2.00a2LM, "24 Jul 2019", an ALPHA 2 release -- has NO
# --check-sex, --impute-sex, or --het flag at all. The originally planned
# one-line `plink2 --check-sex` call from an earlier draft of this script
# does not exist on this build, so this step instead exports the same
# raw ingredients --check-sex would have used internally:
#   - --freq: per-variant allele frequency, for the expected-
#     heterozygosity-under-HWE side of the calculation.
#   - --export A: per-sample additive (0/1/2/NA) genotypes, for the
#     observed-heterozygosity side.
# r_scripts/sex_check_viz.R combines these into the same X inbreeding
# coefficient (F) --check-sex computes: F = 1 - (observed het / expected
# het under HWE), with expected het per variant = 2*p*(1-p). This is the
# standard method-of-moments inbreeding-coefficient formula, not a guess
# at --check-sex's internals.
#
# --split-par hg38 IS confirmed present and spelled correctly in this
# build's own --help text ("'b38'/'hg38' = GRCh38, 2781479/155701383")
# -- unlike --check-sex, this one needed no fix. Run as a separate plink2
# call from the immediately following --chr X restriction (rather than
# combined into one call), since whether --chr filtering sees pre- or
# post-split chromosome codes within a single invocation isn't something
# this build's docs make explicit -- two small, individually-verifiable
# steps sidesteps the ambiguity rather than risking a silent wrong
# answer.
#
# Re-imports chrX from the cohort VCF directly, rather than reusing
# relatedness/cohort_qc.* (step 25): that fileset is --autosome-only --
# chrX kinship needs per-sample sex, which is exactly what this step
# exists to determine, so using it as an input here would be circular.
# Same --set-all-var-ids/--new-id-max-allele-len handling as
# plink_relatedness.sh, for the same reason (real structural indels can
# exceed --set-all-var-ids's built-in allele-length cap). --maf 0.05
# --geno 0.05 (no --mind -- this step needs a call for every sample, not
# to drop poorly-genotyped ones) mirrors plink_relatedness.sh's QC
# rationale, for a stable F estimate.

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load plink/2.001

VCF="vqsr/cohort.pass.normalized.vcf.gz"
OUT_DIR="sex_check"

mkdir -p "${OUT_DIR}"

echo "Importing chrX genotypes from the cohort VCF"

plink2 \
  --threads 8 \
  --vcf "${VCF}" \
  --double-id \
  --max-alleles 2 \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 1000 truncate \
  --chr X \
  --make-pgen \
  --out "${OUT_DIR}/cohort_chrX"

echo "Splitting pseudoautosomal regions"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_chrX" \
  --split-par hg38 \
  --make-pgen \
  --out "${OUT_DIR}/cohort_chrX_split"

echo "Restricting to non-PAR X"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_chrX_split" \
  --chr X \
  --make-pgen \
  --out "${OUT_DIR}/cohort_chrX_nonpar"

echo "Filtering low-quality/rare variants (--maf, --geno) for a stable F estimate"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_chrX_nonpar" \
  --maf 0.05 \
  --geno 0.05 \
  --make-pgen \
  --out "${OUT_DIR}/cohort_chrX_qc"

echo "Computing allele frequencies"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_chrX_qc" \
  --freq \
  --out "${OUT_DIR}/cohort_chrX_qc"

echo "Exporting per-sample additive genotypes"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_chrX_qc" \
  --export A \
  --out "${OUT_DIR}/cohort_chrX_qc"

echo "Done -- run r_scripts/sex_check_viz.R to compute per-sample X inbreeding coefficient (F)"
