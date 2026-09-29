#!/bin/bash

# Gathers jobs/samtools_sex_check.sh's per-sample rows into one table.
# Plain shell, run manually once all 121 array tasks have finished --
# same role as make_cohort_map_genomedbi.sh/make_crosscheck_params.sh/
# make_haplotype_sites_bed.sh elsewhere in this repo.

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

OUT_DIR="sex_check"
OUT_FILE="${OUT_DIR}/cohort_sex_depth_ratios.tsv"

echo -e "sample\tchr1_mapped\tchrX_mapped\tchrY_mapped\tchrX_ratio\tchrY_ratio" > "${OUT_FILE}"
cat "${OUT_DIR}"/*.sex_depth.txt >> "${OUT_FILE}"

n_rows=$(($(wc -l < "${OUT_FILE}") - 1))
echo "Wrote ${n_rows} sample rows to ${OUT_FILE}"
