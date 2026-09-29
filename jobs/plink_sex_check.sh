#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name plink_sex_check
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 16G
#SBATCH --time 2:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A.log
#SBATCH --verbose

# Confirms each WGS donor's biological sex from X-chromosome heterozygosity,
# as a QC checkpoint against the sample_demographics.csv sex column and
# against jobs/samtools_sex_check.sh's independent chrX/chrY depth-ratio
# check. Not a SLURM array -- one set of chrX genotype calls across the
# whole cohort at once, same shape as plink_relatedness.sh.
#
# Males are hemizygous for X outside the pseudoautosomal regions (PAR1/
# PAR2, which are present on both X and Y and so remain diploid in males
# too), so their non-PAR X genotype calls should come back essentially
# all-homozygous; females should show normal autosome-like heterozygosity.
# plink2 --check-sex computes an X inbreeding coefficient (F) per sample
# from this and calls sex from it (F near 1 => male, F near 0 => female,
# an ambiguous band in between flagged for manual review).
#
# IMPORTANT, UNVERIFIED against this cluster's plink/2.001 module: this
# repo's own history (see ancestry_pca.sh) has already hit real cases
# where this specific build -- self-reports internally as "24 Jul 2019"
# -- predates PLINK2 features/syntax that current web docs assume, twice
# needing a fix based on this build's own `plink2 --help <flag>` output
# rather than generic docs. Both --split-par and --check-sex below are
# used with current-docs syntax but have NOT been confirmed against this
# build's own --help text the way ancestry_pca.sh's --pca/--score calls
# eventually were. If either errors out, run `plink2 --help split-par`
# and/or `plink2 --help check-sex` on the cluster and share the output
# before guessing again -- that's what worked for --pca.
#
# Re-imports chrX from the cohort VCF rather than reusing
# relatedness/cohort_qc.* (step 25): that fileset is --autosome-only
# (chrX kinship needs per-sample sex, which is exactly what this step
# exists to determine -- using it as an input here would be circular).
# Same --set-all-var-ids/--new-id-max-allele-len handling as
# plink_relatedness.sh, for the same reason (real structural indels in a
# 121-genome cohort can exceed --set-all-var-ids's built-in allele-length
# cap).
#
# --split-par hg38 relabels the PAR1/PAR2 portion of chrX (and chrY, if
# present) to a separate pseudo-chromosome before the sex check, so those
# diploid-in-both-sexes regions don't dilute the non-PAR X-heterozygosity
# signal --check-sex actually needs. Confirm "hg38" is this build's
# correct keyword for GRCh38 PAR boundaries via --help if this errors --
# some PLINK2 versions expect "b38" instead.
#
# --maf 0.05 --geno 0.05 (no --mind, deliberately -- this step needs a
# sex call for every sample, not to drop poorly-genotyped ones) mirrors
# plink_relatedness.sh's QC rationale: a stable X-heterozygosity/F
# estimate benefits from filtering out rare and poorly-genotyped variants
# first, standard practice in published sex-check protocols (e.g.
# Anderson et al. 2010).

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

echo "Filtering low-quality/rare variants (--maf, --geno) for a stable F estimate"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_chrX_split" \
  --maf 0.05 \
  --geno 0.05 \
  --make-pgen \
  --out "${OUT_DIR}/cohort_chrX_qc"

echo "Running the X-heterozygosity sex check"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_chrX_qc" \
  --check-sex \
  --out "${OUT_DIR}/cohort_sex_check"

echo "Done"
echo "Confirm the actual output filename/extension and column layout --"
echo "PLINK1.9's --check-sex wrote a .sexcheck file with columns FID IID"
echo "PEDSEX SNPSEX STATUS F; plink2's may differ. r_scripts/sex_check_viz.R"
echo "prints whatever it actually finds before assuming column names."
