#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name ancestry_pca
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 32G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A.log
#SBATCH --verbose

# Estimates genetic ancestry for each WGS sample by projecting the cohort
# onto a PCA computed from the 1000 Genomes Phase 3 (hg38) reference
# panel -- intended as a covariate for QTL mapping. Requires
# relatedness/cohort_qc.* (from plink_relatedness.sh) to already exist,
# and the 1000G reference panel already downloaded/decompressed at
# REF_PFILE below -- confirm that path/prefix matches what's actually on
# disk before running.
#
# Projects the cohort onto the reference panel's PCA (--score against
# PLINK2's --pca biallelic-var-wts loadings) rather than merging the two
# datasets with --pmerge: --pmerge was still under active development as
# of 2022, and the plink/2.001 module here self-reports as a build from
# 24 Jul 2019 (confirmed in an earlier job's log on this cluster) --
# solidly before --pmerge existed, so a --pmerge-based approach would
# likely fail outright. Projection also avoids needing to reconcile
# REF/ALT allele coding and strand orientation across a full dataset
# merge, and is widely considered the more rigorous approach for this
# kind of reference-panel ancestry inference anyway. Follows PLINK2's own
# documented recipe: https://www.cog-genomics.org/plink/2.0/score#pca_project
#
# Both sides get re-IDed to a common chrom:pos:ref:alt scheme before
# intersecting -- the cohort VCF and 1000G don't share a variant-ID
# convention. --new-id-max-allele-len 20 missing is safe here (unlike in
# plink_relatedness.sh, which needed "truncate" mode to preserve
# uniqueness for indels) because --snps-only excludes every indel in the
# same call, so no over-length allele ever survives into the final ID'd
# output -- nothing that would hit that length cap is left to collide.
#
# A wider LD-pruning window than plink_relatedness.sh's kinship pruning
# (200 50 0.1) is deliberate, not a typo: PCA is much more sensitive to
# residual LD than KING-robust kinship is, so PCA-oriented pruning
# conventionally uses a wider window.

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load plink/2.001

# CONFIRM this matches the actual downloaded/decompressed 1000G filenames
# under /projects/p31535/boles/plink_references/ before running.
REF_PFILE="/projects/p31535/boles/plink_references/all_hg38"
COHORT_PFILE="relatedness/cohort_qc"
OUT_DIR="ancestry"

mkdir -p "${OUT_DIR}"

echo "Harmonizing reference panel variant IDs"

plink2 \
  --threads 8 \
  --pfile "${REF_PFILE}" \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 20 missing \
  --autosome \
  --max-alleles 2 \
  --snps-only just-acgt \
  --rm-dup exclude-all \
  --make-pgen \
  --out "${OUT_DIR}/ref_reid"

echo "Harmonizing cohort variant IDs"

plink2 \
  --threads 8 \
  --pfile "${COHORT_PFILE}" \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 20 missing \
  --autosome \
  --max-alleles 2 \
  --snps-only just-acgt \
  --rm-dup exclude-all \
  --make-pgen \
  --out "${OUT_DIR}/cohort_reid"

echo "Intersecting to the SNP set shared by both"

# Two passes, one per side: --extract only ever accepts a plain ID list
# (a .snplist), not a .pvar table, so each side's snplist has to be
# written out explicitly before extracting it from the other side. The
# reference is filtered to the cohort's SNPs first, and its own
# resulting snplist -- not the cohort's original one -- is what then
# gets extracted from the cohort, guaranteeing both sides end up with
# the exact same shared set rather than two lists that merely overlap.
plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_reid" \
  --write-snplist \
  --out "${OUT_DIR}/cohort_snps"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/ref_reid" \
  --extract "${OUT_DIR}/cohort_snps.snplist" \
  --make-pgen \
  --out "${OUT_DIR}/ref_shared"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/ref_shared" \
  --write-snplist \
  --out "${OUT_DIR}/ref_shared_snps"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_reid" \
  --extract "${OUT_DIR}/ref_shared_snps.snplist" \
  --make-pgen \
  --out "${OUT_DIR}/cohort_shared"

echo "LD-pruning the reference panel's shared SNP set"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/ref_shared" \
  --maf 0.05 \
  --indep-pairwise 1000 100 0.1 \
  --out "${OUT_DIR}/ref_pruned"

echo "Computing reference panel PCA and allele frequencies"

plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/ref_shared" \
  --extract "${OUT_DIR}/ref_pruned.prune.in" \
  --freq \
  --pca 10 biallelic-var-wts \
  --out "${OUT_DIR}/ref_pca"

echo "ref_pca.eigenvec.allele header (sanity check: ID should be column 2, A1 column 6, PC1 column 7):"
head -1 "${OUT_DIR}/ref_pca.eigenvec.allele"

echo "Projecting cohort samples onto the reference PCA"

# --read-freq uses the REFERENCE panel's allele frequencies (from --freq
# above), not the cohort's own -- variance-standardize needs to
# standardize against the same population the PCA space was built from,
# not re-derive frequencies from the (much smaller, non-representative)
# cohort being projected.
plink2 \
  --threads 8 \
  --pfile "${OUT_DIR}/cohort_shared" \
  --extract "${OUT_DIR}/ref_pruned.prune.in" \
  --read-freq "${OUT_DIR}/ref_pca.afreq" \
  --score "${OUT_DIR}/ref_pca.eigenvec.allele" 2 6 header-read no-mean-imputation variance-standardize \
  --score-col-nums 7-16 \
  --out "${OUT_DIR}/cohort_projected_pca"

echo "Done"
