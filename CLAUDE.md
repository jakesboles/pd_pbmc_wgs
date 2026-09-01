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
- Job scripts submitted from `/projects/b1169/boles/pd_pbmc_wgs`
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
| `bowtie_params_id.txt` | one row per **sample** (unique ID prefix before first `_`) | sample ID only | `bwa_merge.sh`, `gatk_markduplicates.sh`, `samtools_qc.sh`, `gatk_baserecalibrator.sh`, `gatk_haplotype_caller.sh` |
| `bowtie_params_r1.txt` / `_r2.txt` | one row per sample | comma-joined list of that sample's R1/R2 files across all lanes | unused — written for the now-removed bowtie2 path (see "Removed: legacy bowtie2 path" below) |
| `cutadapt_params.txt` | one row per **lane-level fastq pair** | `R1_file,R2_file,replicate_id` | `cutadapt.sh` |
| `bwa_params.txt` | one row per lane-level fastq pair | `R1_file,R2_file,replicate_id,lane,sample` | `bwa.sh` |

Despite the `bowtie_*` naming, `bowtie_params_id.txt` is the de facto
**master sample manifest** used throughout the bwa/GATK production path —
this is a naming artifact from an earlier alignment approach, not a sign
that bowtie2 is actually in use (see "Removed: legacy bowtie2 path" below).

Cohort scale as encoded in current array sizes: **1804** raw fastq files
(902 lane-level R1/R2 pairs) collapsing down to **121** unique samples.

Unlike the large biological outputs (BAMs, VCFs — see `.gitignore`), these
small text manifest/parameter files are committed to the repo as they're
generated, since they're needed to reproduce or re-run any given step.
`bwa_params.txt`, `cutadapt_params.txt`, `chromosomes.txt` (the
chr1-22+chrX list consumed by steps 13-14), and `cohort.sample_map` (the
`sample<TAB>gvcf_path` map built by step 12, see below) are all currently
checked in this way. `bowtie_params_id.txt`/`_r1.txt`/`_r2.txt` are not
checked in — `cohort.sample_map` now serves as the more up-to-date,
121-sample master list for anything that needs it going forward (e.g. the
crosscheck-fingerprinting crosswalk below).

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

21. **`jobs/make_crosscheck_params.sh`** — build the WGS↔scATAC crosswalk
    that `gatk_crosscheckfingerprints.sh` needs. Plain shell, run manually
    (not a SLURM array), same role as `make_cohort_map_genomedbi.sh`.
    Repurposes `cohort.sample_map` as the WGS sample list; for each WGS
    sample it strips the `JSB` prefix to get the Cell Ranger ARC output
    directory name (e.g. `JSB100-1` → `100-1`), then reads the real `SM`
    tag out of that directory's `atac_possorted_bam.bam` header (rather
    than assuming it matches the directory name).
    In: `cohort.sample_map`,
    `/projects/b1042/Gate_Lab/boles/pd_pbmc_multiome/cellranger/<code>/outs/atac_possorted_bam.bam`.
    Out: `crosscheck_sample_map.txt` (`wgs_sample<TAB>atac_sample`, only
    for samples with a matching multiome directory and a readable `SM`
    tag), `crosscheck_atac_bams.txt` (one matched BAM path per line, same
    order/count as the sample map), `crosscheck_missing_atac.txt` (WGS
    samples with no matching Cell Ranger ARC directory).

22. **`jobs/gatk_crosscheckfingerprints.sh`** — `gatk
    CrosscheckFingerprints`, comparing each WGS sample's genotypes in the
    cohort VCF against its matched scATAC BAM's genotype-likelihood
    signal at haplotype-map SNP sites, to confirm donor identity between
    the two datasets. Uses `--INPUT_SAMPLE_MAP` to rename each VCF sample
    to its scATAC `SM` tag for comparison (so the `JSB`-prefix mismatch
    doesn't block matching), `--CROSSCHECK_BY SAMPLE` (the GATK default is
    `READGROUP`, which would compare below the level we want), and
    `--EXIT_CODE_WHEN_MISMATCH 0` so a genotype mismatch — a real possible
    QC finding, e.g. a sample swap — doesn't get treated as a job failure.
    In: `vqsr/cohort.pass.normalized.vcf.gz`, `crosscheck_sample_map.txt`,
    `crosscheck_atac_bams.txt`,
    `/projects/p31535/boles/Homo_sapiens_assembly38.haplotype_database.txt`.
    Out: `crosscheck/cohort_vs_atac.crosscheck_metrics` — per-sample-pair
    `LOD_SCORE` and `RESULT` (e.g. `EXPECTED_MATCH`,
    `EXPECTED_MISMATCH`) to review before trusting any WGS↔multiome sample
    pairing downstream.

### Removed: legacy bowtie2 path

`jobs/bowtie2_build.sh`, `jobs/bowtie2.sh`, and the dev/test
`jobs/gatk_haplotype_caller_test.sh` (which read from the bowtie2 output
directory `bam/` instead of the production `bwa_bam/*.bqsr.bam`) have been
removed from the repo — `workflow.txt` never referenced this path, and the
production aligner has always been `bwa-mem2` (step 6). `bowtie_params_r1.txt`/
`_r2.txt`, which only ever fed `bowtie2.sh`, are correspondingly unused now
(see "Sample/lane bookkeeping" above). Note `multiqc_config.yaml`'s module
list still includes a `bowtie2` section; that's now stale and can be
dropped whenever MultiQC gets its own SLURM script.

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

## Multiome data layout

Cell Ranger ARC output for the matched scATAC/scRNA (multiome) data lives
at `/projects/b1042/Gate_Lab/boles/pd_pbmc_multiome/cellranger` — note this
stays on `b1042`; only the WGS working directory itself moved to `b1169`
(see "Compute environment" above) — one subdirectory per multiome library,
named with the sample code minus the WGS `JSB` prefix (e.g. `100-1` for
WGS sample `JSB100-1`). Each
subdirectory follows the standard `cellranger-arc count` layout; the files
relevant here are `outs/atac_possorted_bam.bam` (+ `.bai`), used for
fingerprinting, and eventually `outs/gex_possorted_bam.bam` and the
`filtered_feature_bc_matrix*`/`atac_fragments.tsv.gz` outputs for the QTL
mapping stage itself.

## Next step

Run `jobs/make_crosscheck_params.sh` then `jobs/gatk_crosscheckfingerprints.sh`
(steps 21-22 above) and review `crosscheck/cohort_vs_atac.crosscheck_metrics`
for any `EXPECTED_MATCH` sample pair that actually comes back as a
mismatch (or vice versa) before trusting any WGS↔multiome sample pairing.
Two things to double check before/while running:

- The haplotype map path
  (`/projects/p31535/boles/Homo_sapiens_assembly38.haplotype_database.txt`)
  — confirm this exact path and filename exist on `p31535` as written.
- `crosscheck_missing_atac.txt` (written by `make_crosscheck_params.sh`) —
  any WGS sample listed there has no matching Cell Ranger ARC directory
  and won't be checked; confirm whether that's expected (e.g. a WGS-only
  sample with no multiome library) or a naming mismatch to fix.

Once identity is confirmed, the actual QTL mapping work (integrating
`vqsr/cohort.pass.normalized.vcf.gz` genotypes with the multiome
scRNA/scATAC data) has no scripts in this repo yet.
