#!/bin/bash
# Ad hoc diagnostic -- run manually, not part of the pipeline.
# Finds which ATAC BAM(s) have a different reference-contig order than
# the rest of the cohort (the cause of the CrosscheckFingerprints
# SequenceListsDifferException).

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load samtools/1.16.1-gcc-10.4.0

OUT="params/crosscheck_bam_dict_signatures.txt"
> "$OUT"

paste params/crosscheck_sample_map.txt params/crosscheck_atac_bams.txt | while IFS=$'\t' read -r wgs_sample atac_sample bam; do
  sig=$(samtools view -H "$bam" \
    | awk -F'\t' '/^@SQ/ { for (i = 1; i <= NF; i++) if ($i ~ /^SN:/) print $i }' \
    | md5sum | cut -d' ' -f1)
  echo -e "${wgs_sample}\t${atac_sample}\t${sig}\t${bam}" >> "$OUT"
done

echo "=== Signature counts (majority = expected/consistent ordering) ==="
cut -f3 "$OUT" | sort | uniq -c | sort -rn

echo
echo "=== Samples NOT matching the majority signature ==="
majority_sig=$(cut -f3 "$OUT" | sort | uniq -c | sort -rn | head -1 | awk '{print $2}')
awk -F'\t' -v maj="$majority_sig" '$3 != maj' "$OUT"