#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name genesis_pcrelate
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 32G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A.log
#SBATCH --verbose

# Prep step for the ancestry-aware relatedness re-analysis via GENESIS's
# PC-Relate -- a more rigorous alternative to plink_relatedness.sh's plain
# KING-robust kinship (step 25), which doesn't account for population
# structure at all. This job only exports genotypes; the actual PC-Relate
# run is r_scripts/genesis_pcrelate.R (run afterward via
# `Rscript r_scripts/genesis_pcrelate.R`, or interactively). Requires
# r_scripts/install_genesis_packages.R to have been run interactively
# first (see that script for why it can't just be part of a batch job),
# and relatedness/cohort_qc.* and relatedness/cohort_pruned.prune.in
# (from plink_relatedness.sh, step 25) to already exist.
#
# PC-Relate needs the pruned genotypes in classic PLINK1 BED/BIM/FAM
# format (SNPRelate::snpgdsBED2GDS, used to build the GDS file GENESIS
# operates on, doesn't read PLINK2's .pgen), so this re-exports exactly
# the same pruned marker set already used for cohort_king.kin0 --
# extracting from cohort_qc.pgen, not the full cohort_raw, so the same
# variant/sample QC already applied for kinship estimation carries over
# here too.
#
# See r_scripts/genesis_pcrelate.R for the actual PC-Relate logic. Note
# this no longer uses GENESIS's PC-AiR or kingToMatrix() at all (an
# earlier version did, and needed several fixes for it -- see git history
# if curious) -- confirmed against the GENESIS source
# (UW-GAC/GENESIS R/pcrelate.R) that pcrelate() only needs a plain
# character vector of sample IDs for `training.set` and any numeric
# matrix with sample-ID rownames for `pcs`, neither requiring a pcair()
# object. Since this repo already has validated ancestry PCs
# (ancestry/cohort_ancestry_pcs_corrected.tsv, from ancestry_pca.sh +
# ancestry_check_scoring.sh/ancestry_viz.R) and hand-picks its own
# training set directly, PC-AiR has nothing left to contribute here. The
# training set itself is chosen by farthest-point sampling directly on
# those ancestry PCs (no kinship data involved), not from
# cohort_king.kin0 -- an earlier version excluded samples by raw KING
# kinship first, but that kinship is exactly what's ancestry-biased here,
# and using it to gate candidacy just starved the training set of the
# minority-ancestry samples that most needed representation.

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load plink/2.001

OUT_DIR="genesis"

mkdir -p "${OUT_DIR}"

echo "Exporting the pruned, QC'd cohort genotypes to PLINK1 BED/BIM/FAM"

plink2 \
  --threads 8 \
  --pfile relatedness/cohort_qc \
  --extract relatedness/cohort_pruned.prune.in \
  --make-bed \
  --out "${OUT_DIR}/cohort_pruned"
