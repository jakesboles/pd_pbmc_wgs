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
  - `plink_references/` — the 1000 Genomes Phase 3 (hg38) PLINK2
    reference panel (genotyped, not sites-only, with population labels),
    used by `ancestry_pca.sh` for ancestry inference. Downloaded manually
    from PLINK2's own resources page, not scripted in this repo.
- Key module/tool versions: `fastqc/0.12.0`, `cutadapt/4.2`, `bwa-mem2
  2.2.1` (built locally at `/projects/p31535/boles/bwa-mem2-2.2.1_x64-linux`,
  not a module), `samtools/1.16.1-gcc-10.4.0`, `gatk/4.4.0.0`,
  `bowtie2/2.5.4`, `R/4.4.0`, `plink/2.001` (self-reports internally as
  build "24 Jul 2019" — notably predates some newer PLINK2 features like
  `--pmerge`; see `ancestry_pca.sh`). The VQSR scripts hardcode a Java 17
  binary path rather than trusting the `gatk` module's bundled JRE.

## Sample/lane bookkeeping

All small parameter/manifest files described below (and referenced
throughout this doc) live under `params/` in the working directory, not
the working directory root — moved there in a repo-wide reorganization;
every job script that reads or writes one of them was updated to the
`params/<file>` path accordingly.

`jobs/make_job_params.R` is the hub that turns the raw `fastq/` directory
listing into every downstream job's parameter file. It is **run manually
(not via SLURM)** whenever the fastq directory changes, and produces:

| File | Grain | Columns | Consumed by |
|---|---|---|---|
| `params/bowtie_params_id.txt` | one row per **sample** (unique ID prefix before first `_`) | sample ID only | `bwa_merge.sh`, `gatk_markduplicates.sh`, `samtools_qc.sh`, `gatk_baserecalibrator.sh`, `gatk_haplotype_caller.sh` |
| `params/bowtie_params_r1.txt` / `_r2.txt` | one row per sample | comma-joined list of that sample's R1/R2 files across all lanes | unused — written for the now-removed bowtie2 path (see "Removed: legacy bowtie2 path" below) |
| `params/cutadapt_params.txt` | one row per **lane-level fastq pair** | `R1_file,R2_file,replicate_id` | `cutadapt.sh` |
| `params/bwa_params.txt` | one row per lane-level fastq pair | `R1_file,R2_file,replicate_id,lane,sample` | `bwa.sh` |

Despite the `bowtie_*` naming, `params/bowtie_params_id.txt` is the de facto
**master sample manifest** used throughout the bwa/GATK production path —
this is a naming artifact from an earlier alignment approach, not a sign
that bowtie2 is actually in use (see "Removed: legacy bowtie2 path" below).

Cohort scale as encoded in current array sizes: **1804** raw fastq files
(902 lane-level R1/R2 pairs) collapsing down to **121** unique samples.

Unlike the large biological outputs (BAMs, VCFs — see `.gitignore`), these
small text manifest/parameter files are committed to the repo as they're
generated, since they're needed to reproduce or re-run any given step.
`params/bwa_params.txt`, `params/cutadapt_params.txt`, `params/chromosomes.txt` (the
chr1-22+chrX list consumed by steps 13-14), and `params/cohort.sample_map` (the
`sample<TAB>gvcf_path` map built by step 12, see below) are all currently
checked in this way. `params/bowtie_params_id.txt`/`_r1.txt`/`_r2.txt` are not
checked in — `params/cohort.sample_map` now serves as the more up-to-date,
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
   In: `params/cutadapt_params.txt` (array 1-902), `fastq/`.
   Out: `trimmed_fastqs/*.fastq.gz`, per-sample log in `cutadapt_logs/`.

4. **`jobs/trimmed_fastqc.sh`** — FastQC on every trimmed fastq file.
   In: `trimmed_fastqs/*.fastq.gz` (array 1-1804).
   Out: `trimmed_fastqc_reports/`.

5. **`jobs/make_job_params.R`** (re-run, or reuse `params/bwa_params.txt` from
   step 2) — builds the BWA-specific parameter file with per-lane sample
   and lane identifiers needed for read-group tagging.

6. **`jobs/bwa.sh`** — align each lane-level trimmed fastq pair to hg38
   with `bwa-mem2 mem`, streamed directly into `samtools sort`. Read group
   is set per lane: `ID=<sample>_<lane>`, `SM=<sample>`,
   `LB=<sample>_lib1`, `PL=ILLUMINA`, `PU=<lane>`.
   Prerequisite (one-time): **`jobs/bwa_build.sh`** — `bwa-mem2 index`,
   `samtools faidx`, `gatk CreateSequenceDictionary` on the reference.
   In: `params/bwa_params.txt` (array 1-902), `trimmed_fastqs/`.
   Out: `bwa_bam/<replicate_id>.sorted.bam` (+ `.bai`) — **one BAM per
   lane**, not yet per sample.

7. **`jobs/bwa_merge.sh`** — merge all lane-level BAMs belonging to one
   sample, then coordinate-sort the merged file.
   In: `params/bowtie_params_id.txt` (sample list), `bwa_bam/<sample>_*.sorted.bam`.
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
    Out: `params/cohort.sample_map`.

13. **`jobs/gatk_genomicsdbiimport.sh`** — combine all samples' GVCFs into
    a per-chromosome GenomicsDB workspace.
    In: `params/cohort.sample_map`, `params/chromosomes.txt` (array 1-23; one chromosome
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
    Repurposes `params/cohort.sample_map` as the WGS sample list; for each WGS
    sample it strips the `JSB` prefix to get the Cell Ranger ARC output
    directory name (e.g. `JSB100-1` → `100-1`), then reads the real `SM`
    tag out of that directory's `atac_possorted_bam.bam` header (rather
    than assuming it matches the directory name).
    In: `params/cohort.sample_map`,
    `/projects/b1042/Gate_Lab/boles/pd_pbmc_multiome/cellranger/<code>/outs/atac_possorted_bam.bam`.
    Out: `params/crosscheck_sample_map.txt` (`wgs_sample<TAB>atac_sample`, only
    for samples with a matching multiome directory and a readable `SM`
    tag), `params/crosscheck_atac_bams.txt` (one matched BAM path per line, same
    order/count as the sample map), `params/crosscheck_missing_atac.txt` (WGS
    samples with no matching Cell Ranger ARC directory).

22. **`jobs/make_haplotype_sites_bed.sh`** — build a BED file of the
    haplotype map's SNP positions. Plain shell, run manually, once. The
    haplotype map is a Picard-format text file: a SAM-style `@HD`/`@SQ`
    header block, then a `#`-prefixed column-header line
    (`#CHROMOSOME  POSITION  NAME  MAJOR_ALLELE  MINOR_ALLELE  MAF  ...`),
    then tab-separated SNP rows; this just pulls `CHROMOSOME`/`POSITION`
    into 0-based BED coordinates.
    In: `/projects/p31535/boles/Homo_sapiens_assembly38.haplotype_database.txt`.
    Out: `params/haplotype_sites.bed`.

23. **`jobs/subset_reorder_atac_bams.sh`** — fixes a reference-contig-order
    mismatch between the ATAC BAMs and the WGS side before
    `CrosscheckFingerprints` can compare them. The Cell Ranger ARC
    reference lists contigs alphabetically (`chr1, chr10, chr11, ...`);
    the WGS cohort VCF and the haplotype map (both built from Broad's
    `Homo_sapiens_assembly38.fasta`) list them numerically
    (`chr1, chr2, chr3, ...`). Same contigs, same lengths, different
    order — but `CrosscheckFingerprints` requires every file it
    fingerprints together to share an identical sequence dictionary, so
    every single WGS/ATAC comparison fails with htsjdk's
    `SequenceListsDifferException` without this step. This affects **all**
    samples, not an isolated one or two — confirmed by hashing all 121
    ATAC BAMs' `@SQ` orderings (all identical to each other) and then
    comparing a representative one directly against the VCF and haplotype
    map (both numeric, differing from the BAMs). A SLURM array, one task
    per `params/crosscheck_sample_map.txt`/`params/crosscheck_atac_bams.txt` line.
    Reordering each full ~15-20GB ATAC BAM via `gatk ReorderSam` would
    mean rewriting/re-sorting the entire file, so each task first subsets
    its BAM down to just the reads overlapping `params/haplotype_sites.bed` (fast
    via the existing `.bai` index — `CrosscheckFingerprints` never looks
    at anything else anyway), *then* reorders that small subset with
    `ReorderSam -SD <WGS .dict>`. Also passes
    `--ALLOW_INCOMPLETE_DICT_CONCORDANCE true`: the subset BAM still
    inherits its full original header, including ALT/unplaced-scaffold
    contigs (e.g. `KI270728.1`) that Broad's `.dict` doesn't declare the
    same way — `ReorderSam` validates *every* header contig against the
    target dictionary regardless of whether any reads actually use it, so
    without this flag it refuses outright on the first unmatched one, even
    though the `-L` subsetting already means no such reads are present.
    In: `params/haplotype_sites.bed`, `params/crosscheck_sample_map.txt`,
    `params/crosscheck_atac_bams.txt`,
    `/projects/p31535/boles/Homo_sapiens_assembly38.dict` (built by
    `bwa_build.sh`, step 6).
    Out: `crosscheck/atac_subset/<wgs_sample>.subset.reordered.bam` (+
    index), one pair per sample.

24. **`jobs/gatk_crosscheckfingerprints.sh`** — `gatk
    CrosscheckFingerprints`, comparing *every* WGS sample's genotypes in
    the cohort VCF against *every* scATAC BAM's genotype-likelihood
    signal at haplotype-map SNP sites — a full all-pairs (121×121)
    comparison, not just the 121 presumed-matched pairs — to definitively
    confirm donor identity between the two datasets and catch any sample
    label swap. **Reworked from an earlier per-pair design** (one SLURM
    array task per presumed WGS/ATAC pair, `--INPUT`/`--SECOND_INPUT`
    each subsetted to just that one pair, `CROSSCHECK_MODE
    CHECK_SAME_SAMPLE`) that could only confirm or deny each presumed
    pairing in isolation: if a sample's true genetic match were actually
    a *different* ATAC BAM than the one presumed, that design would never
    even put the two in the same comparison to notice. One job, not an
    array: `--INPUT` is the full cohort VCF (all 121 samples) with no
    per-sample `SelectVariants` subsetting; `--SECOND_INPUT` is repeated
    once per sample for all 121 *reordered subset* BAMs from step 23
    (`crosscheck/atac_subset/<wgs_sample>.subset.reordered.bam`, built
    from `params/crosscheck_sample_map.txt`), not the raw
    `atac_possorted_bam.bam` — see step 23 for why those need reordering.
    `--INPUT_SAMPLE_MAP params/crosscheck_sample_map.txt` still renames
    each WGS VCF sample to its scATAC `SM` tag so the `JSB`-prefix
    mismatch doesn't block matching. `--CROSSCHECK_MODE CHECK_ALL_OTHERS`
    (not the default `CHECK_SAME_SAMPLE`) is what makes this an all-pairs
    check: it does everything `CHECK_SAME_SAMPLE` does (confirms each
    sample's presumed match) *and* additionally confirms every sample
    does not unexpectedly match any other sample, producing a `RESULT`
    column of `EXPECTED_MATCH`/`EXPECTED_MISMATCH` (presumed pairing
    behaved as expected) or `UNEXPECTED_MATCH`/`UNEXPECTED_MISMATCH`
    (it didn't) for every one of the 121×121 comparisons — an
    `UNEXPECTED_MATCH` between two different presumed samples is exactly
    the sample-swap signature to look for. `--CROSSCHECK_BY SAMPLE` (GATK
    default is `READGROUP`) and `--EXIT_CODE_WHEN_MISMATCH 0` (a genotype
    mismatch, including a real swap, is an expected QC finding to review,
    not a task failure) are unchanged from the old design. Because both
    `INPUT` and `SECOND_INPUT` now carry the full 121-sample cohort on
    each side, the old design's reason for per-sample `INPUT` subsetting
    (avoiding ~120 harmless "sample X is missing from RIGHT group" log
    lines per task, from a 121-vs-1 `INPUT`/`SECOND_INPUT` size mismatch)
    no longer applies — every sample now has a real counterpart on both
    sides. Heavier than the old per-pair tasks (one job computing
    ~14,641 pairwise LOD scores instead of 121 jobs of one each), so
    resources were bumped to 32G/8h from the old 16G/2h.
    In: `vqsr/cohort.pass.normalized.vcf.gz`,
    `params/crosscheck_sample_map.txt`,
    `crosscheck/atac_subset/<wgs_sample>.subset.reordered.bam` (step 23,
    all 121), `/projects/p31535/boles/Homo_sapiens_assembly38.haplotype_database.txt`.
    Out: `crosscheck/cohort_all_pairs.crosscheck_metrics` — one row per
    VCF-sample × BAM-sample comparison (121×121), each with a `LOD_SCORE`
    and `RESULT` — review every `UNEXPECTED_MATCH`/`UNEXPECTED_MISMATCH`
    row before trusting any WGS↔multiome sample pairing downstream; the
    *absence* of any `UNEXPECTED_MATCH` across the full matrix is what
    actually rules out a label swap, which the old per-pair design could
    not do. **`r_scripts/crosscheckfingerprint_scores_viz.R`** reads this
    file directly (`skip = 6` for the metrics header), tabulates
    `RESULT`, and plots a `LEFT_GROUP_VALUE` × `RIGHT_GROUP_VALUE`
    `LOD_SCORE` heatmap (`crosscheck/lod_heatmap.png`) plus a `LOD_SCORE`
    histogram — the heatmap is the fastest way to eyeball an off-diagonal
    hotspot.

24b. **`jobs/make_crosscheck_params_gex.sh`** / **`jobs/subset_reorder_gex_bams.sh`**
    / **`jobs/gatk_crosscheckfingerprints_gex.sh`** — the same all-pairs
    identity check as step 24, applied to each sample's GEX (scRNA)
    BAM (`outs/gex_possorted_bam.bam`, in the same Cell Ranger ARC output
    directory as the ATAC BAM) instead of the ATAC BAM, as an additional,
    independent sanity check on WGS↔multiome donor identity. Mechanically
    identical to steps 21+23+24: `make_crosscheck_params_gex.sh` builds a
    parallel crosswalk (`params/crosscheck_gex_sample_map.txt`,
    `params/crosscheck_gex_bams.txt`, `params/crosscheck_missing_gex.txt`)
    the same way `make_crosscheck_params.sh` does, reading each GEX BAM's
    real RG `SM` tag rather than assuming it matches the ATAC side's (in
    practice expected to match, since both come from the same Cell Ranger
    ARC `--id`, but confirmed fresh rather than assumed);
    `subset_reorder_gex_bams.sh` subsets each GEX BAM to
    `params/haplotype_sites.bed` and reorders it to the WGS reference's
    contig order, same as `subset_reorder_atac_bams.sh` (the same
    alphabetical-vs-numeric mismatch is expected, since one Cell Ranger
    ARC run shares a single reference bundle across its ATAC and GEX
    outputs); `gatk_crosscheckfingerprints_gex.sh` is the same single
    all-pairs `CHECK_ALL_OTHERS` job as step 24's reworked design, just
    pointed at the GEX crosswalk and reordered BAMs.
    **Caveat specific to this comparison:** the haplotype map's SNP sites
    were chosen assuming roughly uniform genome-wide coverage, a
    reasonable assumption for WGS and (open-chromatin-biased but still
    genome-wide) ATAC fragments, but not for RNA-seq — GEX reads only
    cover transcribed, predominantly exonic sequence, and depth at any
    given site tracks that gene's expression level rather than genomic
    position. Expect meaningfully fewer informative (covered) sites per
    sample, and correspondingly weaker (though still decisive for a true
    match, assuming enough covered sites remain) `LOD_SCORE`s than the
    ATAC comparison — an expected property of comparing against RNA-seq,
    not a bug. No RNA-seq-specific preprocessing (e.g. `SplitNCigarReads`)
    is applied, since `CrosscheckFingerprints` does its own pileup-based
    comparison rather than calling variants.
    In: `vqsr/cohort.pass.normalized.vcf.gz`, `params/cohort.sample_map`,
    `params/haplotype_sites.bed`,
    `/projects/p31535/boles/Homo_sapiens_assembly38.dict`,
    `/projects/p31535/boles/Homo_sapiens_assembly38.haplotype_database.txt`,
    `/projects/b1042/Gate_Lab/boles/pd_pbmc_multiome/cellranger/<code>/outs/gex_possorted_bam.bam`.
    Out: `params/crosscheck_gex_sample_map.txt`,
    `params/crosscheck_gex_bams.txt`, `params/crosscheck_missing_gex.txt`,
    `crosscheck/gex_subset/<wgs_sample>.subset.reordered.bam` (+ index),
    `crosscheck/cohort_all_pairs_gex.crosscheck_metrics` — review the
    same way as step 24's ATAC output, keeping the weaker-LOD-score
    caveat above in mind. `r_scripts/crosscheckfingerprint_scores_viz.R`
    can be pointed at this file in place of the ATAC one (swap its
    hardcoded input path, and the heatmap's output filename/axis label)
    to get the same `RESULT` tabulation, LOD heatmap, and histogram for
    the GEX comparison.

25. **`jobs/plink_relatedness.sh`** — a *different* QC axis from steps
    21-24: cryptic relatedness *between WGS subjects themselves*
    (duplicate enrollments, unreported family relationships), not
    identity between the WGS and scATAC datasets. Not a SLURM array —
    one set of pairwise comparisons across the whole cohort at once, same
    shape as `gatk_filter_split.sh`. Imports the cohort VCF into PLINK2
    binary (`.pgen`) format (`--double-id`, `--set-all-var-ids` since the
    VCF's ID column is unannotated — with `--new-id-max-allele-len 1000
    truncate`, since some structural indels' allele strings exceed
    `--set-all-var-ids`'s small built-in length cap, and the "truncate"
    mode keeps IDs unique for the later `--extract` step, unlike the
    "missing" mode PLINK2's own error message suggests, which would give
    every over-length variant the same `.` ID — `--autosome`-only since
    chrX kinship needs per-sample sex info this pipeline doesn't track),
    applies
    variant/sample QC (`--maf 0.05 --geno 0.05 --mind 0.1` — deliberately
    *no* `--hwe`, since HWE deviation is expected at real sites in a
    cohort that may contain relatives, which is exactly what this step is
    checking for), LD-prunes to an approximately independent marker set
    (`--indep-pairwise 200 50 0.1`), then estimates pairwise kinship with
    PLINK2's KING-robust estimator (`--make-king-table`) — robust to
    population stratification, unlike classic IBD/PI_HAT. No
    `--king-table-filter` is set, so the output covers every pairwise
    comparison (~7260 for 121 samples), not just flagged/related ones.
    `--geno` and `--mind` run as two separate, ordered `plink2` calls
    (variants filtered first, then samples) rather than combined into
    one — PLINK2 computes per-sample missingness against whatever variant
    set is currently loaded, so running `--mind` before `--geno` has
    dropped the worst sites let a routine amount of joint-genotyping
    missingness (any sample can lack a confident call at a site private
    to other samples) drag every sample's apparent missingness rate
    over 10%, removing all 121 samples on the first attempt. Filtering
    variants before samples is standard published GWAS QC practice for
    exactly this reason — the two steps aren't commutative.
    In: `vqsr/cohort.pass.normalized.vcf.gz`.
    Out: `relatedness/cohort_raw.*`, `relatedness/cohort_geno.*`,
    `relatedness/cohort_qc.*` (PLINK2 filesets),
    `relatedness/cohort_pruned.prune.in`/`.prune.out`,
    `relatedness/cohort_king.kin0` — the pairwise kinship table to
    review: KING-scale kinship coefficients (~0.5 duplicate/MZ twin,
    ~0.25 1st-degree, ~0.125 2nd-degree, ~0.0625 3rd-degree, halving each
    step out; conventional midpoint cutoffs for calling a category are
    ~0.354/0.177/0.0884/0.0442) — most pairs should cluster near 0, and
    any unexpectedly elevated pair is worth following up before assuming
    cohort subjects are all unrelated.

26. **`jobs/ancestry_pca.sh`** — estimates genetic ancestry per sample,
    for use as a QTL-mapping covariate, by projecting the cohort onto a
    PCA computed from the 1000 Genomes Phase 3 (hg38) reference panel.
    Deliberately does **not** merge the cohort into the reference with
    `--pmerge` (the intuitive approach, and what an early draft of this
    script attempted): `--pmerge` was still under active development as
    of 2022, while the `plink/2.001` module here self-reports internally
    as a build from **24 Jul 2019** — confirmed from an earlier job's own
    log output — solidly predating `--pmerge`. Instead follows PLINK2's
    documented projection recipe
    (https://www.cog-genomics.org/plink/2.0/score#pca_project): compute
    PCA + allele frequencies on the reference panel alone
    (`--pca var-wts --freq` — no explicit PC count, since this old
    PLINK2 build rejects a count alongside a modifier and 10 is its
    default anyway), then project the cohort onto
    those loadings with `--score` using the reference's own allele
    frequencies (`--read-freq`) rather than the cohort's. This sidesteps
    `--pmerge` entirely and avoids reconciling REF/ALT allele coding and
    strand orientation across a full dataset merge — also generally
    considered the more statistically rigorous approach for this kind of
    reference-panel ancestry inference. Both cohort and reference are
    first re-IDed to a shared `chrom:pos:ref:alt` scheme and restricted
    to biallelic SNPs (`--snps-only just-acgt`) before intersecting —
    `--extract` only accepts a plain ID list, not a `.pvar` table, so
    each side's SNP list is written out explicitly and the *reference's*
    resulting list (not the cohort's original one) is what gets
    extracted from the cohort, so both sides end up with the exact same
    shared set. Both re-ID calls also pass `--allow-extra-chr`: the 1000G
    panel's `.pvar` includes ALT/unplaced-scaffold contigs (e.g.
    `chr1_KI270706v1_random`) that PLINK2 refuses to even load by
    default, before `--autosome` gets a chance to drop them — this flag
    only relaxes that load-time check, `--autosome` still restricts the
    actual output the same way. Pruning before PCA uses a wider window
    (`--indep-pairwise 1000 100 0.1`) than `plink_relatedness.sh`'s
    kinship pruning (`200 50 0.1`) — deliberate, since PCA is far more
    sensitive to residual LD than KING-robust kinship is.
    In: `relatedness/cohort_qc.*` (step 25), the 1000 Genomes reference
    panel at `/projects/p31535/boles/plink_references/` (confirm the
    exact filename prefix before running — assumed to be `all_hg38`,
    matching PLINK2's standard resource naming).
    Out: `ancestry/ref_pca.eigenvec` (reference samples' own PC
    coordinates — join with the panel's `.psam` `SuperPop`/`Population`
    columns for known-ancestry labels), `ancestry/cohort_projected_pca.sscore`
    (this cohort's samples projected into that same PC space) — plot PC1
    vs PC2 (or beyond) from both together to see where each WGS sample
    falls relative to labeled reference populations.

26b. **`jobs/ancestry_check_scoring.sh`** / **`r_scripts/ancestry_viz.R`**
    — a scale-correction and validation follow-up to step 26, added
    because PLINK2's `--score`-based projection and its own native `--pca`
    don't produce PCs on the same numeric scale: `--score` output has to
    be rescaled to land in the same units as `--pca`'s native eigenvectors
    before it's meaningful to compare a projected sample directly against
    the reference's own PCA coordinates (e.g. to overlay cohort samples on
    labeled 1000G SuperPop clusters). Rather than trust a theoretical
    scaling formula (`-2/sqrt(eigenvalue)`, computed in `ancestry_viz.R`
    for reference only) blindly, this validates and derives the actual
    correction empirically: `ancestry_check_scoring.sh` re-runs `--score`
    on the reference panel against its *own* PCA loadings (i.e.
    self-projection — the same samples that produced `ref_pca.eigenvec`
    are re-scored via the projection path used for the cohort), and
    `ancestry_viz.R` regresses each PC's native eigenvector value against
    its self-projected score (`lm(eigenvec ~ score)`); the fitted slope is
    the empirical per-PC correction factor, and the fit's R² is checked
    (warns if any PC's R² < 0.999) before trusting it. That correction is
    then applied to both the reference self-projection and the cohort's
    original projection (`ancestry/cohort_projected_pca.sscore` from step
    26), putting cohort samples on the same PC scale as the reference's
    native PCA — validated visually via PC1-vs-PC2 and PC3-vs-PC4 plots
    with the cohort overlaid on labeled 1000G SuperPop groups.
    In: `ancestry/ref_shared`, `ancestry/ref_pruned.prune.in`,
    `ancestry/ref_pca.afreq`, `ancestry/ref_pca.eigenvec.var` (all from
    step 26), `ancestry/ref_pca.eigenvec`, `ancestry/ref_pca.eigenval`,
    `ancestry/cohort_projected_pca.sscore`.
    Out: `ancestry/ref_selfprojected_pca.sscore` (intermediate),
    `ancestry/ancestry_pc1_pc2.png`, `ancestry/ancestry_pc3_pc4.png`
    (validation plots), and **`ancestry/cohort_ancestry_pcs_corrected.tsv`**
    — the cohort's PC1-PC10, corrected onto the reference's native PCA
    scale, used as the ancestry-PC covariate by step 27 below. Originally
    written to the repo root and committed directly (as a deliberate
    small-covariate-table exception in the same spirit as
    `params/cohort.sample_map`); moved into `ancestry/` and is now
    gitignored like the rest of that directory's outputs, consistent
    with every other computed file here rather than a special case.

27. **`jobs/genesis_pcrelate_prep.sh`** / **`r_scripts/genesis_pcrelate.R`**
    — ancestry-aware relatedness re-analysis using Bioconductor GENESIS's
    PC-Relate, a follow-up to step 25's plain KING-robust kinship (which
    doesn't account for population structure at all). `genesis_pcrelate_prep.sh`
    only does the PLINK1 BED/BIM/FAM export PC-Relate's GDS file needs
    (`SNPRelate::snpgdsBED2GDS()` doesn't read PLINK2's `.pgen`); the
    actual analysis is `r_scripts/genesis_pcrelate.R`, run afterward
    separately (`Rscript r_scripts/genesis_pcrelate.R`, or interactively).
    Requires **`r_scripts/install_genesis_packages.R`** to have been run
    once, interactively — `BiocManager::install()` can prompt
    `Update all/some/none? [a/s/n]:`, which hangs forever in a
    non-interactive job.
    **This cohort's ancestry is skewed** (a EUR-majority bulk plus several
    small, distinct minority-ancestry clusters), which causes plain
    KING-robust kinship to systematically overestimate relatedness
    *within* the minority clusters — samples sharing a rare ancestry
    background also share more alleles by descent-from-population, which
    naive KING kinship can't distinguish from true recent relatedness.
    Several increasingly complicated designs were tried and abandoned
    while working through this (GENESIS's PC-AiR run once, then
    iteratively re-run against refined kinship estimates, with an
    automatic or artificially-loosened relatedness threshold to build its
    "unrelated" training set) before landing on a much simpler one — see
    git history if the earlier attempts are ever relevant. **The current
    design skips PC-AiR entirely.** Confirmed directly against the
    GENESIS source (`UW-GAC/GENESIS` `R/pcrelate.R`): `pcrelate()`'s
    `training.set` argument only needs to be a plain character vector of
    sample IDs (checked only for membership in `sample.include`, no
    `pcair()`-derived class required), and its `pcs` argument only needs
    to be any numeric matrix with sample-ID rownames. This repo already
    has both of the things PC-AiR would otherwise be used to produce:
    validated, ancestry-representative PCs
    (`ancestry/cohort_ancestry_pcs_corrected.tsv`, from step 26 +
    step 26b) and a directly hand-picked unrelated training set (below)
    — so PC-AiR (and `kingToMatrix()`, and the bugs both came with —
    plink2's KING column-naming mismatch, a missing `thresh` collapsing
    the whole cohort into one relatedness cluster) have nothing left to
    contribute. It's arguably more robust for this specific cohort too:
    PC-AiR's own PCs are necessarily derived from the cohort's own
    (ancestry-biased) kinship, while the corrected ancestry PCs used here
    come entirely from an external 1000G reference panel that never sees
    this cohort's relatedness or ancestry skew at all.
    One pass, no iteration:
    1. Build the GDS file PC-Relate operates on (same
       `showfile.gds(closeall = TRUE)` handling as before — `gdsfmt`
       tracks open GDS files in an in-process table independent of the
       filesystem, so a prior run's unclosed handle can make
       `snpgdsBED2GDS()` fail with "has been created or opened" even
       after the `.gds` file itself is deleted; this call is a safe
       no-op otherwise).
    2. Load the corrected, validated ancestry PCs, matched/reordered
       against the GDS's own sample IDs with a hard error on any
       mismatch.
    3. Select an ancestry-diverse training set directly from the
       corrected ancestry PCs via **farthest-point (MaxMin) sampling**:
       greedily add whichever remaining sample is farthest, in
       (unit-variance-scaled) PC space, from everyone already selected,
       so the first picks are the most extreme/outlying ancestry points
       and later picks fill in the rest of the spread. **No kinship data
       used at all.** An earlier version of this step instead excluded
       anyone with elevated raw KING kinship (step 25) from candidacy —
       replaced, not kept, after that approach turned out to be
       self-defeating: raw KING is exactly the signal this whole analysis
       exists to correct for ancestry bias, and gating candidacy on it
       meant the stricter the cutoff, the more it starved candidacy of
       the minority-ancestry samples that most needed representation
       (confirmed empirically — loosening that cutoff substantially
       changed the final kinship estimates). Farthest-point sampling
       sidesteps the bias rather than tuning around it: true close
       relatives (parent-child, full sibs) share ~50% of their genome and
       so sit very near each other in ancestry-PC space, meaning
       maximizing spread naturally disfavors picking two of them
       together, without needing to trust the biased raw numbers at all.
       How many samples to keep (`n_training_samples`) is also
       data-driven by default: farthest-point sampling's per-step "gain"
       (how far the newly-added sample was from everyone already
       selected) shrinks monotonically as the training set fills in —
       large gains early (real outliers), small gains later
       (increasingly redundant, interior points) — so the R script finds
       the elbow in that decreasing curve automatically (standard
       max-distance-from-the-chord method) and uses it unless
       `n_training_samples` is set manually after looking at
       `genesis/training_set_selection_curve.png`. Writes
       `genesis/training_set_pc1_pc2.png`/`_pc3_pc4.png` to check the
       resulting selection by eye against
       `ancestry/ancestry_pc1_pc2.png`/`_pc3_pc4.png`; any visible
       ancestry outlier not covered can be raised via
       `n_training_samples` or added directly to
       `manual_force_include_ids`, then the whole script re-run (cheap
       now that it's one pass).
    4. Run `pcrelate()` **once**, with this hand-picked set as
       `training.set` and the corrected ancestry PCs as `pcs`. This is
       the final answer — no second pass, no re-running unless the
       training-set selection itself needs adjusting (step 3).
    In: `relatedness/cohort_qc.*`, `relatedness/cohort_pruned.prune.in`
    (step 25 — for the GDS conversion only; this script no longer reads
    `cohort_king.kin0` at all, see step 3 above),
    `ancestry/cohort_ancestry_pcs_corrected.tsv` (step 26b).
    Out: `genesis/cohort_pruned.{bed,bim,fam}` (intermediate),
    `genesis/cohort.gds`, `genesis/training_set_selection_curve.tsv`/`.png`
    (the full farthest-point ranking and its diversity-gain curve, with
    the auto-suggested elbow marked — the main diagnostic for picking
    `n_training_samples`), `genesis/cohort_manual_training_set.txt` (the
    final hand-picked sample IDs, one per line),
    `genesis/training_set_pc1_pc2.png`/`_pc3_pc4.png` (selection sanity
    check against the full cohort spread), `genesis/cohort_kinship_pcrelate.tsv`
    — pairwise, ancestry-adjusted kinship, categorized with the same
    thresholds as step 25 (~0.354/0.177/0.0884/0.0442) — and
    `genesis/cohort_kinship_pcrelate.png` (kinship-vs-k0 plot). To pick
    `n_pcs_for_adjustment` in the R script (which also sets the
    dimensionality of the farthest-point sampling above), look at
    `ancestry/ref_pca.eigenval` (the reference panel's own PCA
    eigenvalues — the scale the corrected PCs were fit to) for a
    scree-plot elbow, and step 26b's `ancestry/ancestry_pc1_pc2.png`/
    `ancestry_pc3_pc4.png` for how many PCs still visibly separate 1000G
    SuperPop clusters. The R script also logs
    `pcrelate_result$kinBtwn`'s actual column names before referencing
    `kin`, so a wrong assumption there will be immediately visible in the
    job log rather than silently producing an empty/wrong output.

28. **`jobs/plink_sex_check_prep.sh`** / **`jobs/samtools_sex_check.sh`** +
    **`jobs/make_sex_depth_table.sh`** / **`r_scripts/sex_check_viz.R`**
    — confirms each donor's biological sex, as a QC checkpoint against
    the demographics sheet's self-reported sex (catching sample mix-ups)
    and as a sanity check on WGS↔multiome donor linkage. Two independent
    signals, deliberately not requiring any new tool install (no
    `somalier`):
    - **`plink_sex_check_prep.sh`** — X-chromosome heterozygosity. Males
      are hemizygous for X outside the pseudoautosomal regions (PAR1/
      PAR2, diploid in both sexes), so their non-PAR X genotypes should
      come back essentially all-homozygous; females show normal
      heterozygosity, summarized as an X inbreeding coefficient (F): F
      near 1 (near-zero observed heterozygosity) => male, F near 0
      (normal heterozygosity) => female. **A first draft of this script
      called `plink2 --check-sex` directly, which turned out not to
      exist**: confirmed against this cluster's actual `plink2 --help`
      output (the full top-level flag listing, checked directly rather
      than assumed from current web docs) that this build — self-reports
      as PLINK v2.00a2LM, "24 Jul 2019", an **alpha 2** release — has no
      `--check-sex`, `--impute-sex`, or `--het` flag at all. `--split-par
      hg38`, by contrast, *is* confirmed present and correctly spelled in
      this build's own `--help` text (`'b38'/'hg38' = GRCh38,
      2781479/155701383`) — no fix needed there, same as `--score`'s
      modifiers in `ancestry_pca.sh` turned out to already match current
      docs while `--pca` didn't.
      Since there's no built-in command, this job (now prep-only, split
      the same way as `genesis_pcrelate_prep.sh`) instead exports the raw
      ingredients `--check-sex` would have used internally — per-variant
      allele frequencies (`--freq`) and per-sample additive 0/1/2
      genotypes (`--export A`) — and **`r_scripts/sex_check_viz.R`**
      computes F by hand: `F = 1 - (observed heterozygosity / expected
      heterozygosity under HWE)`, expected heterozygosity per variant =
      `2*p*(1-p)` from its allele frequency — the standard
      method-of-moments inbreeding-coefficient estimator, not a guess at
      `--check-sex`'s internals. Not a SLURM array — one set of chrX
      genotype calls across the cohort at once, same shape as
      `plink_relatedness.sh`. Re-imports chrX from the cohort VCF
      directly, rather than reusing `relatedness/cohort_qc.*` (step 25):
      that fileset is `--autosome`-only — chrX kinship needs per-sample
      sex, which is exactly what this step exists to determine, so using
      it as an input here would be circular. Same
      `--set-all-var-ids`/`--new-id-max-allele-len 1000 truncate`
      handling as `plink_relatedness.sh`, for the same reason (real
      structural indels can exceed `--set-all-var-ids`'s built-in
      allele-length cap). `--split-par` and the following `--chr X`
      restriction run as two separate `plink2` calls, not combined into
      one, since whether `--chr` filtering sees pre- or post-split
      chromosome codes within a single invocation isn't something this
      build's docs make explicit. `--maf 0.05 --geno 0.05` (no `--mind`
      — this step needs a call for every sample, not to drop
      poorly-genotyped ones) mirrors `plink_relatedness.sh`'s QC
      rationale, for a stable F estimate. The `.raw`/`.afreq` column
      layout is confirmed by printing them, not assumed, in
      `sex_check_viz.R`, matching the same discipline used for
      `pcrelate_result$kinBtwn` in step 27.
    - **`samtools_sex_check.sh`** — relative read depth on chrX and chrY,
      normalized against chr1 as an autosomal baseline, from
      `samtools idxstats`. Independent of genotype calling entirely, so
      it works directly off each sample's BAM. Males have roughly half
      the chrX depth of females and non-trivial chrY depth; females have
      essentially zero chrY depth outside the PAR. Uses
      `bwa_bam/<sample>.bqsr.bam` (step 10's fully processed,
      recalibrated BAM — the same one HaplotypeCaller calls from). This
      BAM was aligned against the *full* reference genome (including
      chrY), even though this pipeline's variant calling
      (`params/chromosomes.txt`, steps 13-14) never joint-genotypes
      chrY — that restriction only affects which sites get called into
      the cohort VCF, not which reads got mapped into the BAM, so chrY
      read counts are available here even though chrY appears nowhere
      else in this pipeline. A SLURM array (1-121,
      `params/bowtie_params_id.txt`), same convention as
      `gatk_baserecalibrator.sh`, so one bad/missing BAM only fails that
      task. Each task writes one row to its own file;
      **`make_sex_depth_table.sh`** (plain shell, run manually once all
      121 tasks finish, same role as `make_cohort_map_genomedbi.sh`/
      `make_crosscheck_params.sh`) concatenates them into one table.
    **`sex_check_viz.R`** plots the depth-ratio check (expect two
    clusters: females near `chrX_ratio≈1, chrY_ratio≈0`; males near
    `chrX_ratio≈0.5, chrY_ratio>0` — look at the real spread before
    trusting any specific numeric cutoff, a commented-out example call
    is provided as a starting point, not a validated threshold), computes
    and plots the X-heterozygosity F statistic described above (a
    histogram — again, look at the real distribution rather than
    hardcoding a male/female cutoff), and sets up (but doesn't hardcode a
    path for, since none is documented anywhere in this repo) a join
    against whatever sheet tracks each donor's self-reported/clinical sex
    — any mismatch between self-reported and
    genetic sex is a real flag worth resolving before trusting that
    donor's data downstream, not something to silently pick one source
    over the other on.
    In: `vqsr/cohort.pass.normalized.vcf.gz` (`plink_sex_check_prep.sh`),
    `bwa_bam/<sample>.bqsr.bam` (step 10, `samtools_sex_check.sh`),
    `params/bowtie_params_id.txt`.
    Out: `sex_check/cohort_chrX*` (intermediate PLINK2 filesets),
    `sex_check/cohort_chrX_qc.afreq`/`.raw` (per-variant allele
    frequencies and per-sample additive genotypes — the inputs to
    `sex_check_viz.R`'s F computation),
    `sex_check/cohort_sex_check_fstat.tsv` (`IID`, `n_variants`,
    `obs_het`, `exp_het`, `F` — written by `sex_check_viz.R`),
    `sex_check/<sample>.sex_depth.txt` (per-sample, intermediate),
    `sex_check/cohort_sex_depth_ratios.tsv` (`sample`, `chr1_mapped`,
    `chrX_mapped`, `chrY_mapped`, `chrX_ratio`, `chrY_ratio`).

### Removed: legacy bowtie2 path

`jobs/bowtie2_build.sh`, `jobs/bowtie2.sh`, and the dev/test
`jobs/gatk_haplotype_caller_test.sh` (which read from the bowtie2 output
directory `bam/` instead of the production `bwa_bam/*.bqsr.bam`) have been
removed from the repo — `workflow.txt` never referenced this path, and the
production aligner has always been `bwa-mem2` (step 6). `params/bowtie_params_r1.txt`/
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

## WGS↔scATAC identity crosscheck: previously "complete" under a weaker check, now reworked for a definitive all-pairs recheck

Steps 21-23 (`make_crosscheck_params.sh`, `make_haplotype_sites_bed.sh`,
`subset_reorder_atac_bams.sh`) have been run end-to-end on the full
cohort and their outputs (the crosswalk, the haplotype-sites BED, and
all 121 reordered/subset ATAC BAMs) are reusable as-is. Step 24
(`gatk_crosscheckfingerprints.sh`) was originally run as a per-pair SLURM
array (one task per presumed WGS/ATAC pair, `--INPUT`/`--SECOND_INPUT`
each subsetted to just that one sample, `CROSSCHECK_MODE
CHECK_SAME_SAMPLE`); all 121 presumed pairs reported
`RESULT=EXPECTED_MATCH` with `LOD_SCORE` well above the significance
threshold (>>20 for every sample).

**That result only confirmed each presumed pairing looked internally
consistent — it could not rule out a sample label swap.** Each task's
`--INPUT`/`--SECOND_INPUT` only ever contained the one presumed-matching
sample on each side, so even if (say) `JSB100-1`'s true genetic match
were actually the BAM labeled for `JSB100-2`, that comparison had no
opportunity to notice — the real match was never in the comparison at
all. 121 independent one-pair tasks can only answer "does this presumed
pair match," never "does this sample match anyone *else* instead."
Revisited after noticing issues downstream in the multiome analysis that
suggested a possible demographics/sample-identity problem; `step 24`
above has been reworked into a single job that compares the full
121-sample cohort VCF against all 121 reordered ATAC BAMs at once
(`CROSSCHECK_MODE CHECK_ALL_OTHERS`), producing the complete 121×121
LOD/RESULT matrix so an `UNEXPECTED_MATCH` anywhere — the actual
swap signature — would be visible. **Not yet re-run with this design;
the old per-pair result above should not be treated as having ruled out
a swap.** See step 24 above for the full design and rationale, and
"Next step" below for what to do with the output once it's run.

Background on why steps 22-23 exist, for future reference: the first
attempt at `gatk_crosscheckfingerprints.sh` (as a single non-array job
comparing the VCF against all 121 ATAC BAMs at once) crashed with a
`SequenceUtil$SequenceListsDifferException`. The initial read on that
error was that one ATAC BAM had an odd reference build — but hashing all
121 BAMs' `@SQ` orderings showed they're all identical to each other; the
real mismatch is that *all* of them (Cell Ranger ARC's reference) list
contigs alphabetically while the WGS VCF/haplotype map (Broad's
`Homo_sapiens_assembly38.fasta`) list them numerically. That's a
cohort-wide, not per-sample, issue — every comparison would have failed
the same way. Steps 22-23 fix it once for the whole cohort rather than
per sample, and that output is unaffected by the step-24 rework above —
it's reused as-is for the all-pairs comparison.

**Step 24b** mirrors this entire step-24 design against each sample's
GEX (scRNA) BAM instead of the ATAC BAM — a second, independent
all-pairs identity check, added for the same reason as the step-24
rework (a possible demographics/sample-identity problem surfaced in the
multiome analysis) and not yet run. See step 24b above for the full
design and its RNA-seq-coverage caveat (expect fewer informative sites
and weaker `LOD_SCORE`s than the ATAC comparison — not itself a sign of
a problem).

## Next step

**`jobs/gatk_crosscheckfingerprints.sh`** (step 24, reworked) — submit
this job (no array bounds to set any more, just `sbatch`) to get the
definitive all-pairs WGS↔scATAC identity check described above. Review
`crosscheck/cohort_all_pairs.crosscheck_metrics` for any row with
`RESULT` of `UNEXPECTED_MATCH` or `UNEXPECTED_MISMATCH` — the former is
the direct signature of a sample label swap (a WGS sample matching a
*different* sample's ATAC BAM better than its own presumed one); the
diagonal (presumed pairs) should all still come back `EXPECTED_MATCH`
with high `LOD_SCORE` as before, but the point of this rework is the
off-diagonal cells, which the old per-pair design never computed at all.

**`jobs/make_crosscheck_params_gex.sh`** / **`jobs/subset_reorder_gex_bams.sh`**
/ **`jobs/gatk_crosscheckfingerprints_gex.sh`** (step 24b) — run this
same sequence (crosswalk → subset/reorder → single all-pairs job) for
the GEX BAMs. `subset_reorder_gex_bams.sh`'s `--array` bound needs
setting to match `wc -l params/crosscheck_gex_bams.txt` once
`make_crosscheck_params_gex.sh` has run (may differ slightly from the
ATAC crosswalk's 121 if any sample's Cell Ranger ARC directory is
missing one BAM type but not the other). Review
`crosscheck/cohort_all_pairs_gex.crosscheck_metrics` the same way as the
ATAC output; any `UNEXPECTED_MATCH` here independently corroborates (or
contradicts) whatever the ATAC-based all-pairs check above finds, which
is the point of running both — two independent assays agreeing on donor
identity is stronger evidence than either alone, and if they disagree,
that disagreement itself is worth chasing down. This is directly
relevant to the demographics/sample-identity concerns
raised after adding `sample_demographics.csv` — resolve any unexpected
match here before trusting that `r_scripts/sex_check_viz.R`'s sex-based
cross-check (or anything else keyed on sample identity) is comparing the
right donor's WGS and multiome data to each other.

**`jobs/plink_relatedness.sh`** (step 25) has been debugged to a working
state (two real bugs found and fixed along the way: an allele-length
overflow in `--set-all-var-ids`, and a `--geno`/`--mind` filter-ordering
issue that was removing every sample — see the script's own comments and
git history for both). Review `relatedness/cohort_king.kin0`: most
pairwise kinship coefficients should cluster near 0, and any
unexpectedly elevated pair (see the kinship-scale cutoffs in step 25's
description) is worth following up on before assuming all cohort
subjects are unrelated. `r_scripts/relatedness_viz.R` plots `KINSHIP` vs
`IBS0` from this output.

**`jobs/ancestry_pca.sh`** (step 26) — estimates ancestry per sample by
projecting the cohort onto a 1000 Genomes PCA, for use as a QTL-mapping
covariate. The reference panel filename prefix at
`/projects/p31535/boles/plink_references/` has been confirmed to match
`REF_PFILE` in the script (`all_hg38`), and the pipeline runs cleanly
through ID harmonization, `--allow-extra-chr`, the SNP-set intersection,
and LD-pruning. The `--pca` call needed two rounds of fixing against this
cluster's unusually old (24 Jul 2019) PLINK2 build, whose syntax and
output diverge from current PLINK2 web docs: no numeric PC count
alongside a modifier, the modifier itself is `var-wts` (not
`biallelic-var-wts`, which doesn't exist in this build), and it writes
`.eigenvec.var` (not `.eigenvec.allele`) with columns confirmed directly
from this build's own `plink2 --help pca` output — `ID` col 2, `MAJ` col
3 (used as the `--score` allele column, since `--help pca` confirms PC
signs are computed relative to the major allele), PCs in cols 5-14. The
`--score` call's own modifier keywords (`header-read`,
`no-mean-imputation`, `variance-standardize`) were also checked against
this build's `plink2 --help score` output and, unlike `--pca`, are
spelled/behave exactly as current PLINK2 docs describe — no changes
needed there. This step hasn't been run to completion yet; that's the
next thing to confirm.

**`jobs/ancestry_check_scoring.sh`** / **`r_scripts/ancestry_viz.R`** (step
26b) — added after reviewing step 26's raw output: PLINK2's `--score`
projection and its native `--pca` don't share a numeric scale, so a
correction was needed before cohort samples could be meaningfully
compared against the reference's own PCA coordinates. Derives that
correction empirically (self-project the reference panel, regress native
eigenvectors against self-projected scores per PC, use the fitted slope
as the correction factor, sanity-checked via R² > 0.999) rather than
trusting the theoretical formula outright, then applies it to the
cohort's projection and validates visually against labeled 1000G
SuperPop clusters. See step 26b above for details. Produced
`ancestry/cohort_ancestry_pcs_corrected.tsv`, now used as step 27's ancestry-PC
covariate (see below).

**`jobs/genesis_pcrelate_prep.sh`** / **`r_scripts/genesis_pcrelate.R`**
(step 27) — ancestry-aware relatedness re-analysis via GENESIS PC-Relate.
Went through several increasingly complicated designs while chasing a
real problem: this cohort's ancestry skew causes plain KING-robust
kinship to overestimate relatedness within minority ancestry clusters.
First GENESIS's PC-AiR, run once then iteratively re-run against refined
kinship estimates with an automatic or artificially-loosened relatedness
threshold to build its training set; then, after dropping PC-AiR
entirely (confirmed against the GENESIS source that `pcrelate()` needs
neither a `pcair()`-derived training set nor PC-AiR's own PCs), a
still-kinship-gated design that excluded anyone with elevated raw KING
kinship from training-set candidacy before hand-picking the rest via
k-means. That still had the same underlying problem one level down: raw
KING is exactly the signal being corrected for ancestry bias, so gating
candidacy on it directly was self-defeating — confirmed empirically, since
loosening that cutoff substantially changed the final kinship estimates,
meaning the excluded candidates were mattering.
**Current design: no kinship data used in training-set selection at
all.** The training set is now chosen purely by farthest-point (MaxMin)
sampling directly on the validated, corrected ancestry PCs (step 26b) —
greedily picking whichever remaining sample is most different in
ancestry-PC space from everyone already picked. This works as a
relatedness filter too, not just a diversity one: true close relatives
share ~50% of their genome and so sit very close together in ancestry-PC
space, so maximizing spread naturally avoids picking two of them
together, without needing to trust the biased raw KING numbers to make
that call. Training-set *size* is now data-driven as well, via the elbow
in the sampling's diversity-gain curve, rather than a guessed constant.
See step 27 above for the full design.
Not yet run end-to-end with this design — next step is to submit
`genesis_pcrelate_prep.sh`, run `genesis_pcrelate.R`, look at
`genesis/training_set_selection_curve.png` (confirm the auto-suggested
`n_training_samples` looks reasonable, adjusting it by hand if not), then
eyeball `genesis/training_set_pc1_pc2.png`/`_pc3_pc4.png` against
`ancestry/ancestry_pc1_pc2.png`/`_pc3_pc4.png` for any missed ancestry
outlier (adding it to `manual_force_include_ids` and re-running the whole
script if so — cheap now that it's one pass), then review
`genesis/cohort_kinship_pcrelate.tsv` and
`genesis/cohort_kinship_pcrelate.png` against step 25's plain KING
output.

**`jobs/plink_sex_check_prep.sh`** / **`jobs/samtools_sex_check.sh`**
(step 28) — confirms biological sex per donor via two independent
signals (X-heterozygosity from PLINK2, chrX/chrY depth ratio from
`samtools idxstats`), deliberately without adding a new tool
(`somalier`) to the stack. The `samtools_sex_check.sh` array has been
submitted; `plink_sex_check_prep.sh` needed a real fix first. Its
original design called `plink2 --check-sex` directly — the analyst ran
`plink2 --help` on the actual cluster build and confirmed that flag
(along with `--impute-sex` and `--het`) simply doesn't exist on this
build (PLINK v2.00a2LM, "24 Jul 2019", an alpha 2 release). `--split-par
hg38`, by contrast, checked out fine against the same `--help` output —
no change needed there. Since there's no built-in sex-check command on
this build, `plink_sex_check_prep.sh` was reworked (and renamed, prep-
only now, matching `genesis_pcrelate_prep.sh`'s split) to export the raw
ingredients — `--freq` and `--export A` — and `r_scripts/sex_check_viz.R`
computes the X inbreeding coefficient (F) by hand from those, the
standard method-of-moments formula rather than a guess at what
`--check-sex` would have done internally. Not yet run end-to-end with
this design. Next step: submit `plink_sex_check_prep.sh`, run
`make_sex_depth_table.sh` once the `samtools_sex_check.sh` array
finishes, then `sex_check_viz.R` to compute F, review both signals, and
cross-check against the demographics sheet's self-reported sex — any
mismatch is worth resolving before trusting that donor's data in the
QTL-mapping stage.

Once ancestry, relatedness, and sex are all reviewed, the actual QTL
mapping work (integrating `vqsr/cohort.pass.normalized.vcf.gz` genotypes
with the multiome scRNA/scATAC data) has no scripts in this repo yet.
