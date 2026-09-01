#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name subset_reorder_atac
#SBATCH --nodes 1
#SBATCH --array=1-121
#SBATCH --ntasks-per-node 4
#SBATCH --mem 16G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

# Prepares each matched ATAC BAM for gatk_crosscheckfingerprints.sh.
#
# The ATAC BAMs (Cell Ranger ARC's reference) list contigs in
# alphabetical order (chr1, chr10, chr11, ...); the WGS cohort VCF and
# haplotype map (both built from the Broad Homo_sapiens_assembly38.fasta
# reference) list them numerically (chr1, chr2, chr3, ...). Same contigs,
# same lengths, just a different order -- but CrosscheckFingerprints
# requires every file it fingerprints together to share an identical
# sequence dictionary, so every WGS/ATAC comparison would otherwise fail
# with htsjdk's SequenceListsDifferException. This affects every sample,
# not an isolated one -- confirmed by hashing all 121 ATAC BAMs' @SQ
# orderings and finding them all identical to each other, then comparing
# a representative one against the VCF/haplotype map directly.
#
# Rather than reordering each full ~15-20GB ATAC BAM (gatk ReorderSam
# rewrites every record, so it's a full-file, full-resort operation),
# this first subsets each BAM down to just the reads overlapping the
# haplotype map's fingerprinting SNP sites (haplotype_sites.bed, from
# make_haplotype_sites_bed.sh) -- CrosscheckFingerprints never looks at
# anything else anyway -- then reorders that small subset. Requires
# haplotype_sites.bed and crosscheck_sample_map.txt/crosscheck_atac_bams.txt
# (one array task per line/pair) to already exist.
#
# NOTE: --array bounds above must match `wc -l crosscheck_atac_bams.txt`
# -- update both if the crosswalk is regenerated with a different sample
# count.

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load samtools/1.16.1-gcc-10.4.0
module load gatk/4.4.0.0

WGS_DICT="/projects/p31535/boles/Homo_sapiens_assembly38.dict"
SITES_BED="haplotype_sites.bed"
OUT_DIR="crosscheck/atac_subset"

mkdir -p "${OUT_DIR}"

wgs_sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" crosscheck_sample_map.txt | cut -f1)
bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" crosscheck_atac_bams.txt)

echo "${wgs_sample}"
echo "${bam}"

subset_bam="${OUT_DIR}/${wgs_sample}.subset.bam"
reordered_bam="${OUT_DIR}/${wgs_sample}.subset.reordered.bam"

echo "Subsetting to haplotype-map sites"

samtools view -@ 4 -b -L "${SITES_BED}" -o "${subset_bam}" "${bam}"

echo "Reordering contigs to match the WGS reference dictionary"

# ALLOW_INCOMPLETE_DICT_CONCORDANCE: the subset BAM still inherits its
# full original header from the ATAC BAM, including ALT/unplaced-scaffold
# contigs (e.g. KI270728.1) that ReorderSam validates against the new
# dictionary regardless of whether any actual reads use them. Broad's
# Homo_sapiens_assembly38.dict doesn't declare the same ALT/decoy contig
# set as Cell Ranger ARC's reference, so without this flag ReorderSam
# refuses outright on the first unmatched contig it finds -- even though
# the -L subsetting above already restricted the data itself to
# haplotype-map sites on primary chromosomes only, so no reads on those
# contigs are actually present to drop.
gatk ReorderSam \
  -I "${subset_bam}" \
  -O "${reordered_bam}" \
  -SD "${WGS_DICT}" \
  --ALLOW_INCOMPLETE_DICT_CONCORDANCE true

echo "Indexing reordered subset BAM"

samtools index "${reordered_bam}"

rm "${subset_bam}"

echo "Done"
