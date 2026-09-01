#!/bin/bash
# Run manually (not via SLURM), once. Builds a BED file of the haplotype
# map's SNP positions, used by subset_reorder_atac_bams.sh to pull just
# the fingerprinting loci out of each (huge) ATAC BAM before reordering
# it to match the WGS reference's contig order.
#
# The haplotype map is a Picard-format text file: a SAM-style @HD/@SQ
# header block, then a '#'-prefixed column-header line, then tab-
# separated SNP rows (CHROMOSOME, POSITION, NAME, MAJOR_ALLELE,
# MINOR_ALLELE, MAF, ...) -- only the first two columns are needed here.

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

HAPLOTYPE_MAP="/projects/p31535/boles/Homo_sapiens_assembly38.haplotype_database.txt"
OUT="haplotype_sites.bed"

grep -v '^@' "$HAPLOTYPE_MAP" | grep -v '^#' \
  | awk 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' \
  | sort -k1,1 -k2,2n \
  > "$OUT"

echo "Wrote $(wc -l < "$OUT") sites to ${OUT}"
