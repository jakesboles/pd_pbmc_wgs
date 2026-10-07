#!/bin/bash
# Run manually (not via SLURM), on a login/interactive node, once the
# filtered WGS cohort VCF and the Cell Ranger ARC multiome runs are
# available. Not a SLURM array -- same one-time bookkeeping role as
# make_crosscheck_params.sh, mirrored for the GEX (scRNA) side instead of
# ATAC: an independent sanity check on WGS<->multiome donor identity
# using the companion gex_possorted_bam.bam in the same Cell Ranger ARC
# output directory.
#
# Repurposes params/cohort.sample_map (the WGS sample list already used
# to build the GenomicsDB) as the source of truth for which WGS samples
# exist, and looks up each sample's matching Cell Ranger ARC output
# directory by stripping the "JSB" prefix (e.g. JSB100-1 -> 100-1), same
# convention as the ATAC crosswalk.
#
# Reads the real RG SM tag directly out of each gex_possorted_bam.bam
# header, rather than assuming it matches the directory name or the
# ATAC BAM's SM tag (same sample, same Cell Ranger ARC --id, so in
# practice expected to match the ATAC side's SM tag -- but read fresh
# here rather than assumed, same discipline as the ATAC crosswalk).

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load samtools/1.16.1-gcc-10.4.0

CELLRANGER_DIR="/projects/b1042/Gate_Lab/boles/pd_pbmc_multiome/cellranger"

mkdir -p params

SAMPLE_MAP_OUT="params/crosscheck_gex_sample_map.txt"
BAM_LIST_OUT="params/crosscheck_gex_bams.txt"
MISSING_OUT="params/crosscheck_missing_gex.txt"

> "$SAMPLE_MAP_OUT"
> "$BAM_LIST_OUT"
> "$MISSING_OUT"

while IFS=$'\t' read -r wgs_sample gvcf_path; do
  code="${wgs_sample#JSB}"
  bam="${CELLRANGER_DIR}/${code}/outs/gex_possorted_bam.bam"

  if [[ ! -f "$bam" ]]; then
    echo "${wgs_sample}" >> "$MISSING_OUT"
    continue
  fi

  # Pull the real RG SM tag out of the BAM header rather than assuming it
  # matches the directory name or the ATAC BAM's SM tag.
  gex_sample=$(samtools view -H "$bam" \
    | awk -F'\t' '/^@RG/ { for (i = 1; i <= NF; i++) if ($i ~ /^SM:/) { sub("SM:", "", $i); print $i; exit } }')

  if [[ -z "${gex_sample}" ]]; then
    echo "WARNING: no RG SM tag found in ${bam}, skipping ${wgs_sample}" >&2
    continue
  fi

  echo -e "${wgs_sample}\t${gex_sample}" >> "$SAMPLE_MAP_OUT"
  echo "${bam}" >> "$BAM_LIST_OUT"
done < params/cohort.sample_map

echo "Wrote $(wc -l < "$SAMPLE_MAP_OUT") matched WGS/GEX sample pairs to ${SAMPLE_MAP_OUT}"
echo "Wrote $(wc -l < "$MISSING_OUT") WGS samples with no matching Cell Ranger ARC GEX BAM to ${MISSING_OUT}"
