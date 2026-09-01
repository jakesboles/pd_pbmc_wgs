#!/bin/bash
# Run manually (not via SLURM), on a login/interactive node, once both the
# filtered WGS cohort VCF and the Cell Ranger ARC multiome runs are
# available. Not a SLURM array — this is a one-time bookkeeping step, same
# role as make_cohort_map_genomedbi.sh.
#
# Builds the WGS <-> scATAC crosswalk that gatk_crosscheckfingerprints.sh
# needs. Repurposes cohort.sample_map (the WGS sample list already used to
# build the GenomicsDB) as the source of truth for which WGS samples
# exist, and looks up each sample's matching Cell Ranger ARC output
# directory by stripping the "JSB" prefix (e.g. JSB100-1 -> 100-1), per
# the multiome directory naming convention.
#
# Rather than assume the BAM's RG SM tag matches the directory name, this
# reads it directly out of each atac_possorted_bam.bam header, so the
# crosswalk is correct even if Cell Ranger ARC was run with a different
# --id than the directory suggests.

set -euo pipefail

cd /projects/b1042/Gate_Lab/boles/pd_wgs

module load samtools/1.16.1-gcc-10.4.0

CELLRANGER_DIR="/projects/b1042/Gate_Lab/boles/pd_pbmc_mulitome/cellranger"

SAMPLE_MAP_OUT="crosscheck_sample_map.txt"
BAM_LIST_OUT="crosscheck_atac_bams.txt"
MISSING_OUT="crosscheck_missing_atac.txt"

> "$SAMPLE_MAP_OUT"
> "$BAM_LIST_OUT"
> "$MISSING_OUT"

while IFS=$'\t' read -r wgs_sample gvcf_path; do
  code="${wgs_sample#JSB}"
  bam="${CELLRANGER_DIR}/${code}/outs/atac_possorted_bam.bam"

  if [[ ! -f "$bam" ]]; then
    echo "${wgs_sample}" >> "$MISSING_OUT"
    continue
  fi

  # Pull the real RG SM tag out of the BAM header rather than assuming it
  # matches the directory name.
  atac_sample=$(samtools view -H "$bam" \
    | awk -F'\t' '/^@RG/ { for (i = 1; i <= NF; i++) if ($i ~ /^SM:/) { sub("SM:", "", $i); print $i; exit } }')

  if [[ -z "${atac_sample}" ]]; then
    echo "WARNING: no RG SM tag found in ${bam}, skipping ${wgs_sample}" >&2
    continue
  fi

  echo -e "${wgs_sample}\t${atac_sample}" >> "$SAMPLE_MAP_OUT"
  echo "${bam}" >> "$BAM_LIST_OUT"
done < cohort.sample_map

echo "Wrote $(wc -l < "$SAMPLE_MAP_OUT") matched WGS/ATAC sample pairs to ${SAMPLE_MAP_OUT}"
echo "Wrote $(wc -l < "$MISSING_OUT") WGS samples with no matching Cell Ranger ARC directory to ${MISSING_OUT}"
