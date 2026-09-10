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
#   3. Select an ancestry-diverse training set directly from the ancestry
#      PCs via farthest-point (MaxMin) sampling -- greedily add whichever
#      remaining sample is farthest (in PC space) from everyone already
#      selected, so early picks are the most extreme/outlying ancestry
#      points and later picks fill in the rest of the spread. No kinship
#      data used at all: an earlier version of this step instead excluded
#      anyone with elevated raw KING kinship from candidacy, but raw KING
#      is exactly the signal this whole analysis exists to correct for
#      ancestry bias, and using it to gate candidacy turned out to be
#      self-defeating -- the stricter the cutoff, the more it starved
#      candidacy of the minority-ancestry samples that most needed
#      representation (confirmed empirically: loosening that cutoff
#      substantially changed the final kinship estimates, meaning the
#      excluded candidates were mattering). Farthest-point sampling
#      sidesteps this rather than tuning around it: true close relatives
#      (parent-child, full sibs) share ~50% of their genome and so sit
#      very near each other in ancestry-PC space, meaning maximizing
#      spread naturally disfavors picking two of them together, without
#      needing to trust the biased raw numbers to make that call directly.
#      How many samples to keep is picked from the data too: farthest-point
#      sampling's per-step "gain" (how far the newly-added point was from
#      everyone already selected) shrinks monotonically as the training
#      set fills in -- large gains early (real outliers), small gains
#      later (increasingly redundant, interior points) -- so the elbow in
#      that decreasing curve is a principled, visible stopping point,
#      found automatically here (standard max-distance-from-the-chord
#      method) but fully overridable after looking at the plot this
#      script writes out.
#   4. Run pcrelate() once, with this hand-picked set as `training.set`
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

# How many samples to keep in the hand-picked training set. NULL (the
# default) uses the data-driven elbow in the farthest-point sampling
# "gain" curve (see step 3 below and the file header) -- set this to a
# specific integer instead to override that suggestion after looking at
# genesis/training_set_selection_curve.png, e.g. if the elbow looks too
# aggressive/conservative for this cohort.
n_training_samples <- NULL
# Any sample IDs known (from eyeballing the plots this script writes) to
# be ancestry outliers that got left out -- fill this in by hand after a
# first look at genesis/training_set_pc1_pc2.png and
# training_set_pc3_pc4.png, then just re-run the whole script (it's one
# pass now, so re-running is cheap).
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

# ---- Step 3: hand-pick an ancestry-diverse training set (farthest-point sampling) ----
# Scaled to unit variance per PC before computing distances, so PC1
# doesn't dominate purely because it explains more raw variance than
# later PCs -- reasonable here since n_pcs_for_adjustment was already
# chosen (see that variable's comment) to include only PCs that carry
# real ancestry signal, not noise, so weighting them equally is
# weighting real structure equally, not amplifying noise.
pc_cols <- paste0("PC", seq_len(n_pcs_for_adjustment))
scaled_mat <- scale(as.matrix(ancestry_pcs_ordered[, pc_cols]))
rownames(scaled_mat) <- ancestry_pcs_ordered$IID
n_all <- nrow(scaled_mat)

sq_dist_to <- function(mat, point) rowSums(sweep(mat, 2, point)^2)

# Deterministic seed point (farthest from the overall centroid), not a
# random pick -- the whole ranking that follows is then fully
# reproducible with no set.seed() needed.
centroid <- colMeans(scaled_mat)
seed_idx <- which.max(sq_dist_to(scaled_mat, centroid))

order_idx <- integer(n_all)
gain <- rep(NA_real_, n_all)  # gain[i]: min-distance the i-th point added had to the set already selected
order_idx[1] <- seed_idx
min_dist <- sqrt(sq_dist_to(scaled_mat, scaled_mat[seed_idx, ]))
min_dist[seed_idx] <- -Inf

for (i in 2:n_all) {
  nxt <- which.max(min_dist)
  order_idx[i] <- nxt
  gain[i] <- min_dist[nxt]
  d_nxt <- sqrt(sq_dist_to(scaled_mat, scaled_mat[nxt, ]))
  min_dist <- pmin(min_dist, d_nxt)
  min_dist[nxt] <- -Inf
}

selection_curve <- tibble(
  rank = seq_len(n_all),
  IID = rownames(scaled_mat)[order_idx],
  gain = gain
)
write_tsv(selection_curve, file.path(out_dir, "training_set_selection_curve.tsv"))

# Elbow in the (monotonically non-increasing) gain curve: standard
# max-distance-from-the-chord method over ranks 2..n_all (gain[1] is NA
# -- the seed point has no prior set to be "far" from). Purely a
# suggestion -- overridden by setting n_training_samples above if the
# plot below suggests a different cutoff.
y <- gain[-1]
x <- seq_along(y)
xn <- (x - min(x)) / (max(x) - min(x))
yn <- (y - min(y)) / (max(y) - min(y))
x1 <- xn[1]; y1 <- yn[1]; x2 <- xn[length(xn)]; y2 <- yn[length(yn)]
chord_dist <- abs((y2 - y1) * xn - (x2 - x1) * yn + x2 * y1 - y2 * x1) /
  sqrt((y2 - y1)^2 + (x2 - x1)^2)
suggested_n <- which.max(chord_dist) + 1  # +1: y[1] corresponds to a training set of size 2

if (is.null(n_training_samples)) {
  n_training_samples <- suggested_n
  cat("n_training_samples not set -- using elbow-suggested value:", suggested_n, "\n")
} else {
  cat("Using manually-set n_training_samples =", n_training_samples,
      "(elbow suggested", suggested_n, ")\n")
}

p_curve <- ggplot(selection_curve[-1, ], aes(rank, gain)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = suggested_n, linetype = "dashed", color = "red") +
  labs(title = "Farthest-point sampling: diversity gain per added sample",
       subtitle = paste("Dashed line = suggested elbow at rank", suggested_n),
       x = "Samples selected so far", y = "Distance of newly-added sample to selected set") +
  theme_bw()
ggsave(file.path(out_dir, "training_set_selection_curve.png"), p_curve, width = 8, height = 6, dpi = 150)

training_selection_ids <- rownames(scaled_mat)[order_idx[seq_len(n_training_samples)]]

# Visual sanity check against the full cohort's ancestry spread -- saved
# rather than print()'d, since this script may run non-interactively via
# Rscript. Cross-check both plots against ancestry/ancestry_pc1_pc2.png
# and ancestry_pc3_pc4.png: if a visibly distinct ancestry group (e.g. an
# EAS singleton, AFR-leaning points) isn't covered by the red points
# below, either raise n_training_samples or add its sample ID(s) to
# manual_force_include_ids above, then re-run the whole script.
ancestry_pcs_ordered$in_training <- ancestry_pcs_ordered$IID %in% training_selection_ids
p_train_pc12 <- ggplot() +
  geom_point(data = ancestry_pcs_ordered, aes(PC1, PC2), color = "grey70", alpha = 0.5) +
  geom_point(data = filter(ancestry_pcs_ordered, in_training), aes(PC1, PC2), color = "red", size = 3) +
  labs(title = "Hand-picked training set vs. full cohort ancestry spread (PC1 vs PC2)") +
  theme_bw()
ggsave(file.path(out_dir, "training_set_pc1_pc2.png"), p_train_pc12, width = 8, height = 7, dpi = 150)

if (n_pcs_for_adjustment >= 4) {
  p_train_pc34 <- ggplot() +
    geom_point(data = ancestry_pcs_ordered, aes(PC3, PC4), color = "grey70", alpha = 0.5) +
    geom_point(data = filter(ancestry_pcs_ordered, in_training), aes(PC3, PC4), color = "red", size = 3) +
    labs(title = "Hand-picked training set vs. full cohort ancestry spread (PC3 vs PC4)") +
    theme_bw()
  ggsave(file.path(out_dir, "training_set_pc3_pc4.png"), p_train_pc34, width = 8, height = 7, dpi = 150)
}

training_ids <- unique(c(training_selection_ids, manual_force_include_ids))
cat(length(training_ids), "final hand-picked training set:\n")
print(training_ids)
writeLines(training_ids, file.path(out_dir, "cohort_manual_training_set.txt"))

# ---- Step 4: PC-Relate -- one pass, using the hand-picked training set ----
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
