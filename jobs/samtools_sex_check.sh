#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name samtools_sex_check
#SBATCH --nodes 1
#SBATCH --array=1-121
#SBATCH --ntasks-per-node 2
#SBATCH --mem 4G
#SBATCH --time 0:30:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

# Second, independent biological-sex signal, complementing
# jobs/plink_sex_check.sh's X-heterozygosity approach: relative read
# depth on chrX and chrY, normalized against chr1 as an autosomal
# baseline. Males have roughly half the chrX depth of females (one copy
# vs. two) and non-trivial chrY depth; females have essentially zero
# chrY depth outside the pseudoautosomal regions. Unlike the PLINK2
# check, this needs no genotype calls at all -- just mapped-read counts
# from `samtools idxstats`, so it works directly off each sample's BAM
# and doesn't depend on how chrX/chrY variants happen to be called or
# filtered elsewhere in the pipeline.
#
# Uses bwa_bam/<sample>.bqsr.bam (step 10's fully processed, recalibrated
# BAM -- the same one HaplotypeCaller calls variants from) rather than an
# earlier-stage BAM, so depth reflects the same reads used for the rest
# of this pipeline. Critically, this BAM was aligned with bwa-mem2
# against the full reference genome (Homo_sapiens_assembly38.fasta,
# which includes chrY), even though this pipeline's variant calling
# (params/chromosomes.txt, steps 13-14) never joint-genotypes chrY --
# that restriction only affects which sites get called into the cohort
# VCF, not which reads got mapped into the BAM in the first place, so
# chrY read counts are available here even though they aren't anywhere
# downstream in the VCF-based side of this pipeline.
#
# One task per sample (array 1-121, params/bowtie_params_id.txt) so a
# missing/corrupted BAM only fails that one task. Each task writes one
# row to its own file; jobs/make_sex_depth_table.sh (run manually
# afterward) concatenates all 121 into one table.

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load samtools/1.16.1-gcc-10.4.0

PARAMS_FILE="params/bowtie_params_id.txt"
OUT_DIR="sex_check"

mkdir -p "${OUT_DIR}"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PARAMS_FILE}" | cut -f1 -d,)

echo "${sample}"

echo "Computing chrX/chrY depth relative to chr1"

# idxstats columns (no header): ref_name, ref_length, mapped_reads,
# unmapped_reads. chrX_ratio/chrY_ratio = (mapped/length for that contig)
# divided by (mapped/length for chr1) -- a length-normalized relative
# depth, not an absolute coverage figure, which is all a sex call needs.
samtools idxstats "bwa_bam/${sample}.bqsr.bam" | awk -v sample="${sample}" '
  $1 == "chr1" { len1 = $2; m1 = $3 }
  $1 == "chrX" { lenX = $2; mX = $3 }
  $1 == "chrY" { lenY = $2; mY = $3 }
  END {
    chrX_ratio = (m1 > 0 && len1 > 0 && lenX > 0) ? (mX / lenX) / (m1 / len1) : "NA"
    chrY_ratio = (m1 > 0 && len1 > 0 && lenY > 0) ? (mY / lenY) / (m1 / len1) : "NA"
    printf "%s\t%d\t%d\t%d\t%s\t%s\n", sample, m1, mX, mY, chrX_ratio, chrY_ratio
  }
' > "${OUT_DIR}/${sample}.sex_depth.txt"

echo "Done"
