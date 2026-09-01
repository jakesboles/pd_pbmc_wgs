# CLAUDE.md

Guidance for Claude Code (and future readers) working in this repository.

## Project overview

WGS processing pipeline for >100 PBMC samples from people with Parkinson's
disease (PD) and controls. The immediate deliverable is a cohort-level,
hard-filtered, normalized VCF. The larger goal is QTL mapping: integrating
these WGS genotypes with matched single-cell RNA-seq and single-cell
ATAC-seq data generated from the same PBMC samples (multiome).

**Working rule for this repo: do not edit files directly on `main`.** All
changes (including pipeline scripts, docs, configs) go through a branch and
a pull request so the repo owner can review and merge. Treat this repo as
the source of truth for *scripts and workflow documentation only* — raw
data, BAMs, VCFs, and other large outputs are never committed (see
`.gitignore`) and live on the HPC cluster's scratch/project storage.

## Compute environment

- Northwestern **Quest** HPC cluster, SLURM scheduler.
- Job scripts submitted from `/projects/b1042/Gate_Lab/boles/pd_wgs`
  (allocation `b1042`, partition `genomics`, except the VCF-gather step
  which runs under allocation/partition `b1169` for more memory/tmp space).
- Reference genome and resource bundle live on `/projects/p31535/boles`:
  - `Homo_sapiens_assembly38.fasta` (GRCh38, GATK-style contig naming:
    `chr1...chr22, chrX, ...`)
  - `dbsnp_146.hg38.vcf.gz`
  - `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`
  - `hapmap_3.3.hg38.vcf.gz`, `1000G_omni2.5.hg38.vcf.gz`,
    `1000G_phase1.snps.high_confidence.hg38.vcf.gz` (VQSR training/truth
    resources)
- Key module/tool versions: `fastqc/0.12.0`, `cutadapt/4.2`, `bwa-mem2
  2.2.1` (built locally at `/projects/p31535/boles/bwa-mem2-2.2.1_x64-linux`,
  not a module), `samtools/1.16.1-gcc-10.4.0`, `gatk/4.4.0.0`,
  `bowtie2/2.5.4`, `R/4.4.0`. The VQSR scripts hardcode a Java 17 binary
  path rather than trusting the `gatk` module's bundled JRE.

## Sample/lane bookkeeping

`jobs/make_job_params.R` is the hub that turns the raw `fastq/` directory
listing into every downstream job's parameter file. It is **run manually
(not via SLURM)** whenever the fastq directory changes, and produces:

| File | Grain | Columns | Consumed by |
|---|---|---|---|
| `bowtie_params_id.txt` | one row per **sample** (unique ID prefix before first `_`) | sample ID only | `bwa_merge.sh`, `gatk_markduplicates.sh`, `samtools_qc.sh`, `gatk_baserecalibrator.sh`, `gatk_haplotype_caller*.sh`, `bowtie2.sh` |
| `bowtie_params_r1.txt` / `_r2.txt` | one row per sample | comma-joined list of that sample's R1/R2 files across all lanes | `bowtie2.sh` (legacy path, see below) |
| `cutadapt_params.txt` | one row per **lane-level fastq pair** | `R1_file,R2_file,replicate_id` | `cutadapt.sh` |
| `bwa_params.txt` | one row per lane-level fastq pair | `R1_file,R2_file,replicate_id,lane,sample` | `bwa.sh` |

Despite the `bowtie_*` naming, `bowtie_params_id.txt` is the de facto
**master sample manifest** used throughout the bwa/GATK production path —
this is a naming artifact from an earlier alignment approach, not a sign
that bowtie2 is actually in use (see "Known quirks" below).

Cohort scale as encoded in current array sizes: **1804** raw fastq files
(902 lane-level R1/R2 pairs) collapsing down to **121** unique samples.

## Pipeline steps

Each step below: script → what it does → inputs → outputs. Order matches
`workflow.txt`.

1. **`jobs/fastqc.sh`** — FastQC on every raw fastq file.
   In: `fastq/*.fastq.gz` (array 1-1804, one file per task).
   Out: `fastqc_reports/`.

2. **`jobs/make_job_params.R`** — generate all parameter files described
   above from the contents of `fastq/`. Run once per fastq batch, not a
   SLURM array.

3. **`jobs/cutadapt.sh`** — adapter/quality trimming per lane-level R1/R2
   pair. Nextera adapter (`CTGTCTCTTATACACATCT`) trimmed from both reads,
   `--nextseq-trim 20`, `--minimum-length 20`.
   In: `cutadapt_params.txt` (array 1-902), `fastq/`.
   Out: `trimmed_fastqs/*.fastq.gz`, per-sample log in `cutadapt_logs/`.

4. **`jobs/trimmed_fastqc.sh`** — FastQC on every trimmed fastq file.
   In: `trimmed_fastqs/*.fastq.gz` (array 1-1804).
   Out: `trimmed_fastqc_reports/`.

5. **`jobs/make_job_params.R`** (re-run, or reuse `bwa_params.txt` from
   step 2) — builds the BWA-specific parameter file with per-lane sample
   and lane identifiers needed for read-group tagging.

6. **`jobs/bwa.sh`** — align each lane-level trimmed fastq pair to hg38
   with `bwa-mem2 mem`, streamed directly into `samtools sort`. Read group
   is set per lane: `ID=<sample>_<lane>`, `SM=<sample>`,
   `LB=<sample>_lib1`, `PL=ILLUMINA`, `PU=<lane>`.
   Prerequisite (one-time): **`jobs/bwa_build.sh`** — `bwa-mem2 index`,
   `samtools faidx`, `gatk CreateSequenceDictionary` on the reference.
   In: `bwa_params.txt` (array 1-902), `trimmed_fastqs/`.
   Out: `bwa_bam/<replicate_id>.sorted.bam` (+ `.bai`) — **one BAM per
   lane**, not yet per sample.

7. **`jobs/bwa_merge.sh`** — merge all lane-level BAMs belonging to one
   sample, then coordinate-sort the merged file.
   In: `bowtie_params_id.txt` (sample list), `bwa_bam/<sample>_*.sorted.bam`.
   Out: `bwa_bam/<sample>.merged.bam`, `bwa_bam/<sample>.merged.sorted.bam`.
   *Note:* the array in this script is currently `--array=59,99` — only a
   2-sample subset, not the full 121. Confirm with the analyst whether this
   reflects a deliberate small test batch before scaling to the full
   cohort, or whether the array bounds simply haven't been widened yet.

8. **`jobs/gatk_markduplicates.sh`** — `gatk MarkDuplicates` on the merged,
   sorted per-sample BAM.
   In: `bwa_bam/<sample>.merged.sorted.bam` (array currently `59,99`,
   same caveat as step 7).
   Out: `bwa_bam/<sample>.markdup.bam` (+ index),
   `gatk_reports/<sample>.duplicate_metrics.txt`.

9. **`jobs/samtools_qc.sh`** — QC metrics (`flagstat`, `stats`, `idxstats`,
   `coverage`) computed on both the pre-dedup merged BAM and the
   post-markdup BAM, for comparison.
   In: `bwa_bam/<sample>.merged.sorted.bam` and `.markdup.bam` (array
   currently `59,99`).
   Out: `samtools_reports/<sample>.{flagstat,stats,idxstats,coverage}.txt`
   and the `.markdup.*` equivalents. These feed **MultiQC**
   (`multiqc_config.yaml` defines module order: raw FastQC → Cutadapt →
   trimmed FastQC → Bowtie2 → Samtools → GATK), though no `multiqc` SLURM
   script exists yet in `jobs/` — MultiQC appears to be run manually /
   still to be scripted.

10. **`jobs/gatk_baserecalibrator.sh`** — `BaseRecalibrator` +
    `ApplyBQSR` using dbSNP 146 and Mills & 1000G gold-standard indels as
    known-sites.
    In: `bwa_bam/<sample>.markdup.bam` (array **1-121**, full cohort).
    Out: `gatk_reports/<sample>.recal.table`, `bwa_bam/<sample>.bqsr.bam`.

11. **`jobs/gatk_haplotype_caller.sh`** — per-sample germline variant
    calling in GVCF mode.
    In: `bwa_bam/<sample>.bqsr.bam` (array 1-121).
    Out: `haplotype_caller/<sample>.output.g.vcf.gz`.
    *(`jobs/gatk_haplotype_caller_test.sh` is a leftover dev/test variant
    that reads from `bam/<sample>.sorted.bam` — the bowtie2 output
    directory — instead of the production `bwa_bam/*.bqsr.bam`. It is not
    part of the production path; treat it as scratch/debug only.)*

12. **`jobs/make_cohort_map_genomedbi.sh`** — build the GATK sample map
    (`sample<TAB>gvcf_path`) from all per-sample GVCFs. Plain shell, run
    manually.
    In: `haplotype_caller/*.output.g.vcf.gz`.
    Out: `cohort.sample_map`.

13. **`jobs/gatk_genomicsdbiimport.sh`** — combine all samples' GVCFs into
    a per-chromosome GenomicsDB workspace.
    In: `cohort.sample_map`, `chromosomes.txt` (array 1-23; one chromosome
    per task — matches chr1-22 + chrX per the gather step below, so no
    chrY or chrM/MT is called anywhere in this pipeline).
    Out: `genomics_db/chr<N>_db/`.

14. **`jobs/gatk_genotypegvcf.sh`** — joint genotyping per chromosome.
    In: `genomics_db/chr<N>_db` (array 1-23).
    Out: `genotyped_gvcfs/chr<N>.vcf.gz`.

15. **`jobs/gatk_gather_gvcfs.sh`** — concatenate the 23 per-chromosome
    VCFs (chr1-22, chrX) into one cohort VCF. Runs under `b1169` for extra
    memory/tmp headroom.
    In: `genotyped_gvcfs/chr{1..22,X}.vcf.gz`.
    Out: `gathered_genotyped_gvcf/genotyped_cohort.raw.vcf.gz` — **raw,
    unfiltered, joint-genotyped cohort VCF.**

16. **`jobs/gatk_vqsr_recalibrate_snps.sh`** — index the raw cohort VCF,
    then `VariantRecalibrator` in SNP mode using HapMap/Omni/1000G
    (training/truth) and dbSNP (known) resources; annotations `QD MQ
    MQRankSum ReadPosRankSum FS SOR DP`.
    In: `gathered_genotyped_gvcf/genotyped_cohort.raw.vcf.gz`.
    Out: `vqsr/cohort.snps.recal`, `vqsr/cohort.snps.tranches`,
    `vqsr/cohort.snps.plots.R`.

17. **`jobs/gatk_vqsr_apply_snps.sh`** — `ApplyVQSR` in SNP mode,
    truth-sensitivity filter level 99.5.
    In: raw cohort VCF + SNP recal/tranches from step 16.
    Out: `vqsr/cohort.snps_recalibrated.vcf.gz` (SNPs filtered, indels
    still untouched/unfiltered).

18. **`jobs/gatk_vqsr_recalibrate_indels.sh`** — `VariantRecalibrator` in
    INDEL mode on the SNP-recalibrated VCF, using Mills & 1000G (training/
    truth) and dbSNP (known); annotations `QD FS SOR ReadPosRankSum
    MQRankSum DP`.
    In: `vqsr/cohort.snps_recalibrated.vcf.gz`.
    Out: `vqsr/cohort.indels.recal`, `vqsr/cohort.indels.tranches`.

19. **`jobs/gatk_vqsr_apply_indels.sh`** — `ApplyVQSR` in INDEL mode,
    truth-sensitivity filter level 99.0, applied on top of the
    SNP-recalibrated VCF so both filter sets end up on one FILTER column.
    In: `vqsr/cohort.snps_recalibrated.vcf.gz` + indel recal/tranches.
    Out: `vqsr/cohort.recalibrated.vcf.gz` — **fully VQSR-filtered cohort
    VCF (SNP + indel), sites not yet dropped, just flagged.**

20. **`jobs/gatk_filter_split.sh`** — `SelectVariants --exclude-filtered`
    to drop everything that isn't `PASS`, then
    `LeftAlignAndTrimVariants --split-multi-allelics` to normalize
    indels/left-align and split multiallelic records into biallelic ones.
    In: `vqsr/cohort.recalibrated.vcf.gz`.
    Out: `vqsr/cohort.pass.vcf.gz` → **`vqsr/cohort.pass.normalized.vcf.gz`**
    — this is the **final, analysis-ready cohort VCF**: PASS-only,
    biallelic, left-aligned, hg38, GATK contig naming.

### Legacy/unused: bowtie2 path

`jobs/bowtie2_build.sh` and `jobs/bowtie2.sh` build a Bowtie2 index and
align trimmed reads with Bowtie2, writing to `bam/` and
`bowtie2_reports/`. **`workflow.txt` never references this path** — the
production aligner is `bwa-mem2` (step 6). These two scripts, plus the
`bam/`-reading `gatk_haplotype_caller_test.sh`, appear to be an earlier
alignment approach that was superseded by BWA but left in the repo. The
`multiqc_config.yaml` module list still includes a `bowtie2` section,
consistent with this being a genuine (if abandoned) earlier pass rather
than dead code from a typo. Don't build on the bowtie2 outputs for
anything downstream; treat `bwa_bam/*` as the only production alignments.

## Current status (as of last review)

The final output of the WGS-only pipeline is:

```
vqsr/cohort.pass.normalized.vcf.gz
```

a joint-genotyped, VQSR-filtered (PASS only), normalized/biallelic cohort
VCF across chr1-22 and chrX for the full sample set.

Open item to confirm with the analyst: steps 7-9 (`bwa_merge.sh`,
`gatk_markduplicates.sh`, `samtools_qc.sh`) currently have SLURM arrays
restricted to samples `59,99` rather than the full `1-121` used by steps
10-11 onward — worth checking whether the full cohort has actually been
merged/dedup'd/QC'd (and the array bounds just weren't widened in the
committed script) before treating the final VCF as covering all samples.

## Next step

**GATK `CrosscheckFingerprints`**, comparing the filtered WGS cohort VCF
(`vqsr/cohort.pass.normalized.vcf.gz`) against the possorted BAM from the
matched scATAC-seq (multiome) data, per sample. This checks that the donor
identity encoded in each WGS sample's genotypes matches the genotype
signal recoverable from that donor's scATAC reads — a prerequisite for
confidently pairing WGS genotypes to single-cell RNA/ATAC data before QTL
mapping. This will need a haplotype map (e.g. GATK's
`hg38_v0_Homo_sapiens_assembly38.haplotype_database.txt` or equivalent)
and a sample-to-BAM crosswalk between WGS sample IDs and scATAC library
IDs, neither of which exists in this repo yet.
