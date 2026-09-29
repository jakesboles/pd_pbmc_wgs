# Reviews the two independent biological-sex QC signals: PLINK2's
# X-heterozygosity check (jobs/plink_sex_check.sh) and the chrX/chrY
# read-depth-ratio check (jobs/samtools_sex_check.sh +
# jobs/make_sex_depth_table.sh). Run manually/interactively, like
# relatedness_viz.R and ancestry_viz.R -- not a SLURM job.

library(tidyverse)
library(ggplot2)

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

# ---- X-heterozygosity check (PLINK2) ----
# Filename/columns confirmed by print, not assumed -- see
# jobs/plink_sex_check.sh's header comment on why this is unverified
# against this cluster's plink2 build.
plink_fn <- "sex_check/cohort_sex_check.sexcheck"
if (file.exists(plink_fn)) {
  plink_sex <- read_table(plink_fn, show_col_types = FALSE)
  cat(plink_fn, "columns:", paste(names(plink_sex), collapse = ", "), "\n")
  print(plink_sex)
} else {
  cat(plink_fn, "not found -- check jobs/plink_sex_check.sh's own log for",
      "the actual output filename plink2 wrote, and update this path.\n")
  plink_sex <- NULL
}

# ---- Cross-check the two methods against each other ----
# Only meaningful once plink_fn's real ID column name is confirmed above
# (IID in PLINK1.9-style output, but not assumed here) -- update
# plink_id_col below to match before trusting this join.
if (!is.null(plink_sex)) {
  plink_id_col <- "IID"
  combined <- depth %>%
    left_join(plink_sex, by = c("sample" = plink_id_col))
  print(combined)
}

# ---- Cross-check against self-reported/clinical sex ----
# No demographics file path is hardcoded here -- this repo doesn't
# document one. Join `depth` (and/or `plink_sex`) against whatever sheet
# tracks self-reported sex by sample ID, and treat any mismatch as a real
# flag (sample mix-up somewhere upstream, or a genuine biological edge
# case) worth resolving before trusting that donor's data downstream --
# not something to silently prefer one source over the other on.
