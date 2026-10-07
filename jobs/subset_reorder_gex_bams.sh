#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name subset_reorder_gex
#SBATCH --nodes 1
#SBATCH --array=1-121
#SBATCH --ntasks-per-node 4
#SBATCH --mem 16G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

# Prepares each matched GEX (scRNA) BAM for
# gatk_crosscheckfingerprints_gex.sh -- mirrors subset_reorder_atac_bams.sh
# exactly, applied to gex_possorted_bam.bam instead of
# atac_possorted_bam.bam, as an independent additional sanity check on
# WGS<->multiome donor identity.
#
# The GEX BAMs come from the same Cell Ranger ARC reference bundle as the
# ATAC BAMs (one shared --reference for an `cellranger-arc count` run
# covers both the ATAC and GEX outputs), so the same alphabetical-vs-
# numeric contig-order mismatch against the WGS cohort VCF/haplotype map
# (both built from Broad's Homo_sapiens_assembly38.fasta) is expected
# here too -- not re-confirmed per-sample, on the assumption that a
# single Cell Ranger ARC run uses one reference for both BAMs it
# produces; if CrosscheckFingerprints still throws
# SequenceListsDifferException after this step for the GEX side, that
# assumption needs rechecking.
#
# Same subset-then-reorder approach as the ATAC version, for the same
# reason: reordering a full GEX BAM in place (gatk ReorderSam rewrites
# every record) is unnecessary when CrosscheckFingerprints only ever
# looks at reads overlapping the haplotype map's fingerprinting SNP
# sites (params/haplotype_sites.bed, from make_haplotype_sites_bed.sh,
# reused as-is -- it's reference-coordinate-based, not assay-specific).
# Requires params/haplotype_sites.bed and
# params/crosscheck_gex_sample_map.txt/params/crosscheck_gex_bams.txt
# (one array task per line/pair, from make_crosscheck_params_gex.sh) to
# already exist.
#
# NOTE: --array bounds above must match `wc -l params/crosscheck_gex_bams.txt`
# -- update both if the crosswalk is regenerated with a different sample
# count. This may differ from the ATAC crosswalk's count if any sample's
# Cell Ranger ARC directory is missing a GEX BAM that its ATAC BAM has,
# or vice versa.

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load samtools/1.16.1-gcc-10.4.0
module load gatk/4.4.0.0

WGS_DICT="/projects/p31535/boles/Homo_sapiens_assembly38.dict"
SITES_BED="params/haplotype_sites.bed"
OUT_DIR="crosscheck/gex_subset"

mkdir -p "${OUT_DIR}"

wgs_sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" params/crosscheck_gex_sample_map.txt | cut -f1)
bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" params/crosscheck_gex_bams.txt)

echo "${wgs_sample}"
echo "${bam}"

subset_bam="${OUT_DIR}/${wgs_sample}.subset.bam"
reordered_bam="${OUT_DIR}/${wgs_sample}.subset.reordered.bam"

echo "Subsetting to haplotype-map sites"

samtools view -@ 4 -b -L "${SITES_BED}" -o "${subset_bam}" "${bam}"

echo "Reordering contigs to match the WGS reference dictionary"

# ALLOW_INCOMPLETE_DICT_CONCORDANCE: same rationale as
# subset_reorder_atac_bams.sh -- the subset BAM still inherits its full
# original header, including ALT/unplaced-scaffold contigs Broad's .dict
# doesn't declare the same way, which ReorderSam would otherwise refuse
# on outright even though no reads on those contigs survive the -L
# subsetting above.
gatk ReorderSam \
  -I "${subset_bam}" \
  -O "${reordered_bam}" \
  -SD "${WGS_DICT}" \
  --ALLOW_INCOMPLETE_DICT_CONCORDANCE true

echo "Indexing reordered subset BAM"

samtools index "${reordered_bam}"

rm "${subset_bam}"

echo "Done"
