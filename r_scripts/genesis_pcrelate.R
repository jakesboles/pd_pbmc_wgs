# Ancestry-aware relatedness via GENESIS PC-Relate, using a directly
# hand-picked, ancestry-diverse training set -- one pass, no PC-AiR.
#
# This replaces an earlier, considerably more complicated version of this
# script that ran pcair()/kingToMatrix() in multiple iterating passes to
# work around this cohort's ancestry skew (a EUR-majority bulk plus
# several small, distinct minority-ancestry clusters) biasing PC-AiR's
# automatic training-set partition. That complexity turned out to be
# unnecessary: confirmed directly against the GENESIS source
# (UW-GAC/GENESIS R/pcrelate.R) that pcrelate()'s `training.set` argument
# is just a plain character vector of sample IDs (checked only for
# `training.set %in% sample.include`, no class requirement) and its `pcs`
# argument is just any numeric matrix with sample-ID rownames (checked
# only via `is.matrix(pcs)` and `!is.null(rownames(pcs))`) -- neither
# needs to come from a pcair() object at all. Since this repo already has
# both of the things PC-AiR would otherwise be used to produce --
# validated, ancestry-representative PCs (ancestry/
# cohort_ancestry_pcs_corrected.tsv, from ancestry_pca.sh +
# ancestry_check_scoring.sh/ancestry_viz.R) and, below, a directly
# hand-picked unrelated training set -- there's nothing left for pcair()
# to contribute here, and dropping it also drops kingToMatrix() and the
# bugs that came with it (the plink2 column-naming mismatch, the missing
# `thresh` collapsing the whole cohort into one cluster). It's arguably
# more robust for this specific cohort too: PC-AiR's own PCs are
# necessarily derived from the cohort's own kinship, which is exactly
# what ancestry skew was biasing here, while the corrected ancestry PCs
# used below come entirely from an external 1000G reference panel that
# never sees this cohort's relatedness or ancestry skew at all.
#
# The plan, one pass, no iteration:
#   1. Build the GDS file PC-Relate operates on (unchanged from before).
#   2. Load the corrected, validated ancestry PCs (also unchanged).
#   3. Flag confidently-related pairs directly from
#      relatedness/cohort_king.kin0's raw KING kinship (step 25) at a
#      conservative ~3rd-degree cutoff, and exclude them from training-set
#      candidacy. Raw KING is the kinship signal this whole project
#      exists to correct for ancestry bias, but that bias inflates
#      modest/near-zero kinship among same-ancestry samples -- a true
#      close relative (parent-child, siblings, ~0.25 kinship) still
#      clears a conservative cutoff by a wide margin, so this one-shot
#      filter is a reasonable, honestly-imperfect heuristic: it may
#      occasionally exclude a genuinely unrelated pair from a small,
#      homogeneous ancestry cluster, which just costs a few training-set
#      candidates in that cluster, not a wrong final kinship value for
#      that pair (its actual estimate still comes from the same
#      ancestry-adjusted pcrelate() run as everyone else).
#   4. K-means cluster the remaining candidates across the corrected
#      ancestry PCs and pick representative samples from every cluster,
#      guaranteeing ancestry coverage a purely automatic partition could
#      miss or under-represent.
#   5. Run pcrelate() once, with this hand-picked set as `training.set`
#      and the corrected ancestry PCs as `pcs`. This is the final
#      answer -- no second pass, no re-running.

library(GENESIS)
library(GWASTools)
library(SNPRelate)
library(gdsfmt)
library(dplyr)
library(readr)
library(ggplot2)

setwd("/projects/b1169/boles/pd_pbmc_wgs")

out_dir <- "genesis"
dir.create(out_dir, showWarnings = FALSE)

bed_prefix <- file.path(out_dir, "cohort_pruned")
gds_fn <- file.path(out_dir, "cohort.gds")
out_fn <- file.path(out_dir, "cohort_kinship_pcrelate.tsv")

# How many of the corrected ancestry PCs to hand to PC-Relate for
# ancestry adjustment. Pick this by looking at:
#   1. ancestry/ref_pca.eigenval -- the 1000G reference panel's own PCA
#      eigenvalues (the scale these corrected PCs were fit to). Look for
#      the "elbow" where added PCs stop explaining much more variance.
#   2. ancestry/ancestry_pc1_pc2.png and ancestry/ancestry_pc3_pc4.png
#      (from r_scripts/ancestry_viz.R) -- how many PCs still visibly
#      separate distinct 1000G SuperPop clusters, with this cohort's
#      samples overlaid.
n_pcs_for_adjustment <- 4

# Conservative ~3rd-degree cutoff (2^(-9/2) =~ 0.0442, matching the
# categories used throughout this repo) for flagging confidently-related
# pairs to exclude from training-set candidacy -- see the file header for
# why a conservative cutoff on raw (ancestry-biased) KING kinship is still
# a reasonable one-shot filter here.
related_thresh <- 2^(-9/2)

# Cluster count for k-means over the candidate pool, deliberately more
# than the number of visually distinct groups in
# ancestry/ancestry_pc1_pc2.png / ancestry_pc3_pc4.png (about 4-5 here),
# so small/outlier ancestry groups are more likely to land in their own
# cluster instead of being absorbed into the EUR-majority cluster.
n_training_clusters <- 8
# How many representative samples to keep per cluster.
n_per_cluster <- 3
# Any sample IDs known (from eyeballing the plots this script writes) to
# be ancestry outliers that k-means nonetheless failed to select -- fill
# this in by hand after a first look at genesis/training_set_pc1_pc2.png
# and training_set_pc3_pc4.png, then just re-run the whole script (it's
# one pass now, so re-running is cheap). Losing a whole ancestry branch
# to an unlucky clustering draw defeats the purpose of doing this by hand.
manual_force_include_ids <- character(0)

# ---- Step 1: convert the pruned, QC'd cohort genotypes to GDS format ----
# bed_prefix.bed/.bim/.fam is written by jobs/genesis_pcrelate_prep.sh via
# `plink2 --pfile relatedness/cohort_qc --extract
# relatedness/cohort_pruned.prune.in --make-bed` -- the same pruned
# marker set already used to compute cohort_king.kin0.
#
# gdsfmt tracks open GDS files by path in an internal, in-process table for
# as long as the R session lives -- that table is separate from the
# filesystem, so a prior run in the same session (or an interactively
# re-sourced script) that created/opened cohort.gds and didn't reach
# close(gds_reader) at the bottom (e.g. it errored out first) leaves the
# path marked "open" even after the .gds file itself is deleted from disk.
# createfn.gds() then refuses to (re)create it with "has been created or
# opened", and deleting the file has no effect on that in-memory table --
# only closing the handle, or ending the R process, clears it.
# showfile.gds(closeall = TRUE) force-releases anything gdsfmt is tracking
# for this session regardless of how it was orphaned; it's a no-op if
# nothing is open, so this is safe to run unconditionally on every
# invocation, not just after a prior failure.
showfile.gds(closeall = TRUE)
if (file.exists(gds_fn)) file.remove(gds_fn)

snpgdsBED2GDS(
  bed.fn = paste0(bed_prefix, ".bed"),
  bim.fn = paste0(bed_prefix, ".bim"),
  fam.fn = paste0(bed_prefix, ".fam"),
  out.gdsfn = gds_fn
)

gds_reader <- GdsGenotypeReader(filename = gds_fn)
genoData <- GenotypeData(gds_reader)
sample_ids <- getScanID(genoData)
cat("Loaded", length(sample_ids), "samples from", gds_fn, "\n")

# ---- Step 2: load the corrected, reference-projected ancestry PCs ----
# Written by r_scripts/ancestry_viz.R to ancestry/. Sample IDs (IID) are
# matched and reordered against this GDS's own sample_ids -- not just
# assumed to line up -- and any mismatch fails loudly here rather than
# silently misaligning genotypes and PCs inside pcrelate().
ancestry_pcs_fn <- "ancestry/cohort_ancestry_pcs_corrected.tsv"
ancestry_pcs_raw <- read_tsv(ancestry_pcs_fn, show_col_types = FALSE)
cat(ancestry_pcs_fn, "columns:", paste(names(ancestry_pcs_raw), collapse = ", "), "\n")

missing_ids <- setdiff(sample_ids, ancestry_pcs_raw$IID)
if (length(missing_ids) > 0) {
  stop(
    length(missing_ids), " cohort sample(s) from ", gds_fn,
    " are missing from ", ancestry_pcs_fn, ": ",
    paste(head(missing_ids, 10), collapse = ", "),
    if (length(missing_ids) > 10) ", ..." else ""
  )
}

ancestry_pcs_ordered <- ancestry_pcs_raw[match(sample_ids, ancestry_pcs_raw$IID), ]
ancestry_pcs_mat <- as.matrix(select(ancestry_pcs_ordered, starts_with("PC")))
rownames(ancestry_pcs_mat) <- ancestry_pcs_ordered$IID
ancestry_pcs_mat <- ancestry_pcs_mat[, seq_len(n_pcs_for_adjustment), drop = FALSE]

# ---- Step 3: flag confidently-related pairs directly from raw KING kinship ----
# Same header this pipeline's cohort_king.kin0 always has (confirmed from
# the job log, not generic docs): #FID1 ID1 FID2 ID2 NSNP HETHET IBS0
# KINSHIP. read_table() (whitespace-flexible), not read_tsv(), matching
# r_scripts/relatedness_viz.R.
king_raw <- read_table("relatedness/cohort_king.kin0", show_col_types = FALSE)
cat("cohort_king.kin0 columns:", paste(names(king_raw), collapse = ", "), "\n")

flagged_ids <- king_raw %>%
  filter(KINSHIP > related_thresh) %>%
  select(ID1, ID2) %>%
  unlist() %>%
  unique()
cat(length(flagged_ids), "samples excluded from training-set candidacy",
    "(confidently related to someone, raw KINSHIP >", related_thresh, "):\n")
print(flagged_ids)

# ---- Step 4: hand-pick an ancestry-diverse training set ----
candidates <- ancestry_pcs_ordered %>%
  filter(!IID %in% flagged_ids)
cat(nrow(candidates), "candidates remain for training-set selection\n")

pc_cols <- paste0("PC", seq_len(n_pcs_for_adjustment))
scaled_mat <- scale(as.matrix(candidates[, pc_cols]))
rownames(scaled_mat) <- candidates$IID

set.seed(42)
km <- kmeans(scaled_mat, centers = n_training_clusters, nstart = 25)
candidates$cluster <- km$cluster
# Distance to each sample's OWN cluster centroid, in the same scaled space
# k-means actually clustered in (across all n_pcs_for_adjustment PCs, not
# just PC1/PC2) -- using km$centers directly avoids any mismatch between
# the clustering distances and a separately-recomputed mean.
candidates$dist_to_centroid <- sqrt(rowSums((scaled_mat - km$centers[km$cluster, ])^2))

training_selection <- candidates %>%
  group_by(cluster) %>%
  slice_min(order_by = dist_to_centroid, n = n_per_cluster) %>%
  ungroup()
cat(nrow(training_selection), "samples selected across", n_training_clusters, "clusters\n")
print(table(training_selection$cluster))

# Visual sanity check against the full candidate spread -- saved rather
# than print()'d, since this script may run non-interactively via
# Rscript. Cross-check both plots against ancestry/ancestry_pc1_pc2.png
# and ancestry_pc3_pc4.png: if a visibly distinct ancestry group (e.g. an
# EAS singleton, AFR-leaning points) isn't covered by the red points
# below, add its sample ID(s) to manual_force_include_ids above and
# re-run the whole script.
p_train_pc12 <- ggplot() +
  geom_point(data = candidates, aes(PC1, PC2), color = "grey70", alpha = 0.5) +
  geom_point(data = training_selection, aes(PC1, PC2), color = "red", size = 3) +
  labs(title = "Hand-picked training set vs. candidate ancestry spread (PC1 vs PC2)") +
  theme_bw()
ggsave(file.path(out_dir, "training_set_pc1_pc2.png"), p_train_pc12, width = 8, height = 7, dpi = 150)

if (n_pcs_for_adjustment >= 4) {
  p_train_pc34 <- ggplot() +
    geom_point(data = candidates, aes(PC3, PC4), color = "grey70", alpha = 0.5) +
    geom_point(data = training_selection, aes(PC3, PC4), color = "red", size = 3) +
    labs(title = "Hand-picked training set vs. candidate ancestry spread (PC3 vs PC4)") +
    theme_bw()
  ggsave(file.path(out_dir, "training_set_pc3_pc4.png"), p_train_pc34, width = 8, height = 7, dpi = 150)
}

training_ids <- unique(c(training_selection$IID, manual_force_include_ids))
cat(length(training_ids), "final hand-picked training set:\n")
print(training_ids)
writeLines(training_ids, file.path(out_dir, "cohort_manual_training_set.txt"))

# ---- Step 5: PC-Relate -- one pass, using the hand-picked training set ----
# Confirmed against the GENESIS source (UW-GAC/GENESIS R/pcrelate.R):
# `training.set` only needs to be a character vector of sample IDs
# present in `sample.include` (no pcair()-derived class required), and
# `pcs` only needs to be a numeric matrix with sample-ID rownames -- so
# training_ids and ancestry_pcs_mat can be handed to pcrelate() directly.
genoIter <- GenotypeBlockIterator(genoData)

pcrelate_result <- pcrelate(
  gdsobj = genoIter,
  pcs = ancestry_pcs_mat,
  training.set = training_ids,
  sample.include = sample_ids,
  BPPARAM = BiocParallel::SerialParam()
)
cat("pcrelate_result$kinBtwn columns:",
    paste(names(pcrelate_result$kinBtwn), collapse = ", "), "\n")

kin_adjusted <- pcrelate_result$kinBtwn %>%
  mutate(category = case_when(
    kin > 0.354  ~ "Duplicate/MZ twin",
    kin > 0.177  ~ "1st-degree",
    kin > 0.0884 ~ "2nd-degree",
    kin > 0.0442 ~ "3rd-degree",
    TRUE         ~ "Unrelated"
  ))
print(table(kin_adjusted$category))

p_kinship <- ggplot(kin_adjusted, aes(k0, kin)) +
  geom_hline(yintercept = 2^(-seq(3, 9, 2) / 2),
             linetype = "dashed",
             color = "grey") +
  geom_point(alpha = 0.5) +
  theme_bw()
ggsave(file.path(out_dir, "cohort_kinship_pcrelate.png"), p_kinship, width = 8, height = 7, dpi = 150)

write_tsv(kin_adjusted, out_fn)
cat("Wrote", nrow(kin_adjusted), "pairwise estimates to", out_fn, "\n")

close(gds_reader)
