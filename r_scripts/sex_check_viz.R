# Reviews the two independent biological-sex QC signals: PLINK2's
# X-heterozygosity check (jobs/plink_sex_check_prep.sh) and the chrX/chrY
# read-depth-ratio check (jobs/samtools_sex_check.sh +
# jobs/make_sex_depth_table.sh). Run manually/interactively, like
# relatedness_viz.R and ancestry_viz.R -- not a SLURM job.

library(tidyverse)
library(ggplot2)
library(ggbeeswarm)

setwd("/projects/b1169/boles/pd_pbmc_wgs")

# ---- Depth-ratio check (samtools) ----
depth <- read_tsv("sex_check/cohort_sex_depth_ratios.tsv", show_col_types = FALSE)
cat("cohort_sex_depth_ratios.tsv columns:", paste(names(depth), collapse = ", "), "\n")

# Expect two clusters: females near (chrX_ratio ~1, chrY_ratio ~0), males
# near (chrX_ratio ~0.5, chrY_ratio > 0) -- look at the real spread here
# before trusting any specific numeric cutoff; the commented-out call
# below is a starting point, not a validated threshold.
ggplot(depth, aes(chrX_ratio, chrY_ratio, label = sample)) +
  geom_point() +
  labs(title = "chrX/chrY depth relative to chr1, per sample") +
  theme_linedraw()

# depth <- depth %>%
#   mutate(depth_sex_call = case_when(
#     chrY_ratio > 0.1 ~ "Male",
#     chrY_ratio <= 0.1 ~ "Female",
#     TRUE ~ NA_character_
#   ))

# ---- X-heterozygosity check (PLINK2, computed by hand) ----
# This cluster's plink2 build (self-reports "24 Jul 2019", v2.00a2LM --
# an alpha 2 release) has no --check-sex, --impute-sex, or --het flag at
# all, confirmed directly against its own `plink2 --help` output rather
# than assumed from generic docs (see jobs/plink_sex_check_prep.sh's
# header comment). That job instead exports the same raw ingredients
# --check-sex would use internally -- per-variant allele frequencies
# (cohort_chrX_qc.afreq) and per-sample additive genotypes
# (cohort_chrX_qc.raw) -- and this reproduces its X inbreeding
# coefficient by hand:
#   F_i = 1 - (observed heterozygosity for sample i) /
#             (expected heterozygosity for sample i under HWE)
# with expected heterozygosity per variant = 2*p*(1-p) from its allele
# frequency -- the standard method-of-moments inbreeding-coefficient
# estimator, not a guess at --check-sex's internals. F near 1 (near-zero
# observed heterozygosity) => male; F near 0 (normal heterozygosity) =>
# female.
freq_fn <- "sex_check/cohort_chrX_qc.afreq"
raw_fn <- "sex_check/cohort_chrX_qc.raw"

if (file.exists(freq_fn) && file.exists(raw_fn)) {
  freq <- read_tsv(freq_fn, show_col_types = FALSE)
  cat(freq_fn, "columns:", paste(names(freq), collapse = ", "), "\n")

  # .raw files are space-delimited with a literal "NA" for missing calls
  # -- read_table's default na = "NA" handles this without extra options.
  raw <- read_table(raw_fn, show_col_types = FALSE)
  cat(raw_fn, "columns (first 10):",
      paste(head(names(raw), 10), collapse = ", "), "...\n")

  geno_cols <- setdiff(names(raw), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
  # Column names are "<variant ID>_<counted allele>" -- IDs themselves
  # only ever contain colons (--set-all-var-ids '@:#:$r:$a' in
  # plink_sex_check_prep.sh), so splitting off everything after the LAST
  # underscore recovers the ID cleanly, with no ambiguity.
  variant_ids <- sub("_[^_]+$", "", geno_cols)

  # "ID" is confirmed by the print above, not assumed -- update
  # freq_id_col if this build's --freq output names it differently.
  freq_id_col <- "ID"
  p <- freq$ALT_FREQS[match(variant_ids, freq[[freq_id_col]])]
  exp_het_per_variant <- 2 * p * (1 - p)

  geno_mat <- as.matrix(raw[, geno_cols])
  het_mat <- geno_mat == 1  # TRUE where heterozygous; NA propagates for missing calls
  not_na <- !is.na(geno_mat)

  obs_het <- rowSums(het_mat, na.rm = TRUE)
  exp_het <- as.numeric(not_na %*% exp_het_per_variant)

  fstat <- tibble(
    IID = raw$IID,
    n_variants = rowSums(not_na),
    obs_het = obs_het,
    exp_het = exp_het,
    F = 1 - obs_het / exp_het
  )
  write_tsv(fstat, "sex_check/cohort_sex_check_fstat.tsv")
  print(fstat)

  ggplot(fstat, aes(F)) +
    geom_histogram(bins = 40) +
    labs(title = "X inbreeding coefficient (F) per sample",
         subtitle = "F near 1 = male, F near 0 = female") +
    theme_linedraw()
} else {
  cat(freq_fn, "and/or", raw_fn, "not found -- run",
      "jobs/plink_sex_check_prep.sh first.\n")
  fstat <- NULL
}

# ---- Cross-check the two methods against each other ----
if (!is.null(fstat)) {
  combined <- depth %>%
    left_join(fstat, by = c("sample" = "IID"))
  print(combined)
}

# ---- Cross-check against self-reported/clinical sex ----

demographics <- read.csv("sample_demographics.csv")

df <- demographics %>% 
  mutate(sample = paste0("JSB", code)) %>%
  left_join(combined,
            by = "sample")

df %>% 
  ggplot(aes(x = sex,
             y = `F`)) + 
  geom_quasirandom() + 
  theme_linedraw()

df %>% 
  arrange(sex) %>%
  ggplot(aes(x = chrX_ratio,
             y = chrY_ratio)) + 
  geom_point(aes(color = sex),
                   size = 3) + 
  labs(y = "chrY coverage",
       x = "chrX coverage") +
  theme_linedraw()

df %>% 
  filter(sex == "male") %>% 
  arrange(F)
