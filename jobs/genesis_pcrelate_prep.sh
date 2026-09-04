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

# Prep step for the ancestry-adjusted relatedness re-analysis via GENESIS's
# PC-AiR/PC-Relate pipeline -- a more rigorous alternative to
# plink_relatedness.sh's plain KING-robust kinship (step 25), which
# doesn't account for population structure at all. This job only exports
# genotypes; the actual PC-AiR/PC-Relate run is r_scripts/genesis_pcrelate.R
# (run afterward via `Rscript r_scripts/genesis_pcrelate.R`, or
# interactively). Requires r_scripts/install_genesis_packages.R to have
# been run interactively first (see that script for why it can't just be
# part of a batch job), and relatedness/cohort_qc.* and
# relatedness/cohort_pruned.prune.in (from plink_relatedness.sh, step 25)
# to already exist.
#
# PC-AiR needs the pruned genotypes in classic PLINK1 BED/BIM/FAM format
# (SNPRelate::snpgdsBED2GDS, used to build the GDS file GENESIS operates
# on, doesn't read PLINK2's .pgen), so this re-exports exactly the same
# pruned marker set already used for cohort_king.kin0 -- extracting from
# cohort_qc.pgen, not the full cohort_raw, so the same variant/sample QC
# already applied for kinship estimation carries over here too.
#
# See r_scripts/genesis_pcrelate.R for the actual PC-AiR/PC-Relate logic,
# including why PC-Relate's ancestry PCs come from
# cohort_ancestry_pcs_corrected.tsv (jobs/ancestry_check_scoring.sh +
# r_scripts/ancestry_viz.R's corrected 1000G-projection), not PC-AiR's own
# PCs -- and a real bug fixed in an earlier draft of this workflow:
# GENESIS's kingToMatrix() does not understand plink2's --make-king-table
# column names despite what some tutorials assume -- confirmed against the
# GENESIS source rather than guessed.

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
