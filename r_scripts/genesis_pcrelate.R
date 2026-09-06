# Ancestry-adjusted relatedness via GENESIS's PC-AiR / PC-Relate pipeline.
# Run by jobs/genesis_pcrelate_prep.sh (BED/BIM/FAM export step) then this
# script (Rscript r_scripts/genesis_pcrelate.R, or interactively). Requires
# the packages installed by r_scripts/install_genesis_packages.R (run that
# first, interactively).
#
# Note on "ancestry" here: PC-Relate's `pcs` argument below is
# cohort_ancestry_pcs_corrected.tsv -- the ancestry-adjusted, 1000G-
# reference-projected cohort PCs produced by jobs/ancestry_check_scoring.sh
# and r_scripts/ancestry_viz.R (a follow-up correction/validation of
# ancestry_pca.sh's projection, confirmed to separate the cohort cleanly
# along known 1000G SuperPop groups) -- NOT pcair()'s own PCs, even though
# pcair() is still run below (see step 3). This is a deliberate design
# choice, not a shortcut:
#   - PC-AiR's usual selling point is that its PCs aren't confounded by
#     cohort-internal relatedness, unlike a plain PCA run directly on the
#     cohort (where family/duplicate clusters can visibly bias the axes).
#   - The reference-projected PCs used here have that same property for a
#     different reason: their loadings come entirely from the external
#     1000 Genomes reference panel (ancestry_pca.sh's --pca step never
#     sees this cohort at all) -- cohort samples, related or not, are only
#     ever scored/projected onto that fixed external space afterward, so
#     cohort relatedness cannot bias what defines each PC axis.
#   - Combined with ancestry_viz.R's empirical scale correction (matching
#     plink2 --score's projected units to the reference's own native --pca
#     eigenvector units, validated by R^2 > 0.999 self-projection fits)
#     and its confirmed SuperPop separation, these are a defensible
#     substitute for PC-AiR's PCs here -- and arguably preferable, since it
#     keeps kinship estimation on the same ancestry-PC definition likely to
#     be used elsewhere as a QTL-mapping covariate, rather than
#     introducing a second, differently-derived PC basis just for this
#     step.
# pcair() is still run to get a KING-based unrelated "training set" for
# pcrelate() (see step 3) -- that's a different use of PC-AiR than
# supplying PCs, and still needed regardless of PC source.
#
# Design update: this cohort's ancestry is skewed enough (a EUR-majority
# bulk plus several small, distinct minority-ancestry clusters -- see
# ancestry/ancestry_pc1_pc2.png) that plain KING-robust kinship
# systematically overestimates relatedness WITHIN the minority clusters:
# samples sharing a rare ancestry background also share more alleles by
# descent-from-population, which naive KING kinship can't distinguish from
# true recent relatedness. That inflated apparent relatedness caused
# pcair()'s automatic kin.thresh/div.thresh partition (step 3) to exclude
# most of those samples from its "unrelated" training set, leaving it both
# very small and lopsided toward the EUR-majority cluster -- exactly the
# opposite of what PC-Relate needs, and the reason for iterating below:
#   - Step 4/4.5: a first PC-Relate pass, using step 3's (undersized,
#     biased) automatic partition, to get a first-cut, ancestry-corrected
#     kinship estimate (kin_mat_pass1). PC-Relate's ancestry correction is
#     applied to every pair via the `pcs` regression, not just the
#     training-set pairs, so this first pass's kinship numbers are already
#     considerably better-corrected than raw KING, even with a flawed
#     training set.
#   - Step 5: a second, diagnostic-only PC-Relate pass using kin_mat_pass1
#     (rather than raw KING) as pcair()'s kinship input, with GENESIS's own
#     default thresholds (not loosened) -- used only to flag confidently-
#     related pairs for step 6, not as the final answer.
#   - Step 6: instead of further loosening pcair()'s automatic thresholds
#     to force a bigger training set (an earlier attempt at this cluster
#     did exactly that, artificially raising kin.thresh -- a hack that
#     just admits more true near-relatives into "unrelated" rather than
#     fixing the actual problem), manually build an ancestry-diverse
#     unrelated training set: exclude anyone confidently related (step 5),
#     then k-means cluster the remaining candidates across the corrected
#     ancestry PCs and pick representative samples from every cluster --
#     explicitly guaranteeing coverage of small ancestry groups that an
#     automatic threshold-based partition, by chance or by bias, could
#     miss or under-represent.
#   - Step 7/8: pcair()'s own `unrel.set` argument (confirmed against the
#     GENESIS source, UW-GAC/GENESIS R/pcairPartition.R) takes exactly
#     this: it forces the named samples into the unrelated set on top of
#     whatever its automatic kin.thresh/div.thresh partition already
#     found, and -- critically -- never re-flags two samples in unrel.set
#     as "related to each other" even if their pairwise kinship looks
#     elevated, so an ancestry-inflated pair within the hand-picked set
#     isn't second-guessed back out. That's exactly what's needed here;
#     GENESIS's own default thresholds are used again (no more loosening),
#     since the manual set now does the job the threshold hack was trying
#     to hack around. A final PC-Relate pass over this training set
#     produces the kinship estimates actually written out.
#
# PC-AiR needs a preliminary kinship/divergence estimate to find its
# unrelated training set -- that's what plink_relatedness.sh's
# cohort_king.kin0 (already computed) is for. GENESIS::kingToMatrix()
# expects KING-software-style column names (ID1, ID2, Kinship) and does
# NOT recognize plink2's --make-king-table column names directly --
# despite what some tutorials assume, it does not autodetect or support
# the plink2 format (confirmed against the GENESIS source,
# UW-GAC/GENESIS R/makeSparseMatrix.R: it does a strict intersect()
# against literal "ID1"/"ID2"/<estimator> column names). cohort_king.kin0's
# actual header, confirmed by running the job (not guessed from generic
# plink2 docs): `#FID1 ID1 FID2 ID2 NSNP HETHET IBS0 KINSHIP` -- so ID1/ID2
# already match what kingToMatrix wants as-is; only KINSHIP needs renaming
# to Kinship. (An earlier draft of this script assumed IID1/IID2 column
# names, matching plink2's --king-table-format taglist default, and tried
# to rename those -- that assumption was wrong for this build/invocation's
# actual output and would have errored with "can't rename columns that
# don't exist"; fixed here against the real header instead.)

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
king_renamed_fn <- file.path(out_dir, "cohort_king_renamed.kin0")
out_fn <- file.path(out_dir, "cohort_kinship_pcrelate.tsv")

# How many of cohort_ancestry_pcs_corrected.tsv's PCs to hand to PC-Relate
# for ancestry adjustment. The GENESIS vignette's own example uses 2 (for
# PC-AiR PCs), but that's not a default to trust blindly here -- pick this
# by looking at:
#   1. ancestry/ref_pca.eigenval -- the 1000G reference panel's own PCA
#      eigenvalues (the scale these corrected PCs were fit to). Look for
#      the "elbow" where added PCs stop explaining much more variance --
#      PCs past that point are mostly noise, not ancestry structure.
#   2. ancestry/ancestry_pc1_pc2.png and ancestry/ancestry_pc3_pc4.png
#      (written by r_scripts/ancestry_viz.R) -- how many PCs still
#      visibly separate distinct 1000G SuperPop clusters, with this
#      cohort's samples overlaid.
# genesis/cohort_pcair_varprop.txt (from pcair(), step 3 below) is a
# secondary, diagnostic-only cross-check -- it describes PC-AiR's own
# PCs, not the corrected PCs actually used below, but broad agreement
# between the two is a reasonable sanity check that both are picking up
# the same real structure.
n_pcs_for_adjustment <- 4

# Manual training-set selection (step 6). related_thresh is deliberately
# the same conservative ~3rd-degree cutoff used everywhere else in this
# repo (2^(-9/2) =~ 0.0442, matching plink_relatedness.sh's categories) --
# high enough that a merely ancestry-inflated pair is unlikely to cross it
# by chance, so excluding anyone above it from training-set candidacy
# should mostly remove real relatives, not ancestry-skew noise.
related_thresh <- 2^(-9/2)

# Cluster count for k-means over the candidate pool, deliberately more
# than the number of visually distinct groups in
# ancestry/ancestry_pc1_pc2.png / ancestry_pc3_pc4.png (about 4-5 here),
# so that small/outlier ancestry groups are more likely to land in their
# own cluster instead of being absorbed into the EUR-majority cluster.
n_training_clusters <- 8
# How many representative samples to keep per cluster.
n_per_cluster <- 3
# Any sample IDs known (from eyeballing the plots step 6 writes out) to be
# ancestry outliers that k-means nonetheless failed to select -- fill this
# in by hand after a first look at genesis/training_set_pc1_pc2.png and
# training_set_pc3_pc4.png. Losing a whole ancestry branch to an unlucky
# clustering draw defeats the purpose of doing this selection manually.
manual_force_include_ids <- character(0)

# ---- Step 1: convert the pruned, QC'd cohort genotypes to GDS format ----
# bed_prefix.bed/.bim/.fam is written by jobs/genesis_pcrelate_prep.sh via
# `plink2 --pfile relatedness/cohort_qc --extract
# relatedness/cohort_pruned.prune.in --make-bed` -- the same pruned
# marker set already used to compute cohort_king.kin0, so the preliminary
# KING estimate and this re-analysis are on consistent footing.
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

# ---- Step 2: load the existing KING kinship as the preliminary estimate ----
# Confirmed header for this pipeline's cohort_king.kin0 (via plink2
# --make-king-table on the version installed here):
#   #FID1  ID1  FID2  ID2  NSNP  HETHET  IBS0  KINSHIP
# ID1/ID2 already match what kingToMatrix() expects as-is -- only KINSHIP
# needs renaming to Kinship. (FID1/FID2/NSNP/HETHET/IBS0 are left alone;
# kingToMatrix() only reads the columns it needs via intersect() against
# its expected names, so extra columns are harmless.) read_table()
# (whitespace-flexible), not read_tsv(), matching how
# r_scripts/relatedness_viz.R already successfully reads this exact file.
king_raw <- read_table("relatedness/cohort_king.kin0", show_col_types = FALSE)
cat("cohort_king.kin0 columns:", paste(names(king_raw), collapse = ", "), "\n")

king_renamed <- king_raw %>%
  rename(Kinship = KINSHIP)
write_tsv(king_renamed, king_renamed_fn)

# thresh is passed explicitly here -- confirmed against the GENESIS source
# (UW-GAC/GENESIS R/makeSparseMatrix.R): kingToMatrix()'s default is
# thresh = NULL, and with NULL its internal clustering step (used to build
# the sparse block matrix) draws a "relatedness" edge between two samples
# whenever their kinship value is simply != 0 -- not some meaningful
# cutoff. Since plink_relatedness.sh deliberately left --king-table-filter
# unset, cohort_king.kin0 has all ~7260 pairs, including near-zero noise
# values that are nonzero but not remotely "related" -- with thresh=NULL
# every one of those still counts as an edge, collapsing the whole cohort
# into one connected cluster ("121 relatives in 1 clusters; largest
# cluster = 121", "0 samples with no relatives") despite step 25's own
# finding that most pairs cluster near 0 kinship. 2^(-11/2) (~0.0221) is
# GENESIS's own convention for this threshold -- it's the default used by
# kingToMatrix()'s snpgdsIBDClass method, and matches pcair()'s own
# kin.thresh/div.thresh defaults -- so it's used explicitly here too,
# rather than leaving it to the NULL default.
king_mat <- kingToMatrix(
  king_renamed_fn,
  estimator = "Kinship",
  sample.include = sample_ids,
  thresh = 2^(-11/2)
)

# ---- Step 3: PC-AiR -- only used here for its unrelated training set ----
# Uses the same KING matrix for both kinship AND divergence, per GENESIS
# convention: KING-robust kinship already encodes ancestry divergence in
# its negative values. pcair_result$vectors (PC-AiR's own PCs) are written
# out below purely as a diagnostic cross-check against
# cohort_ancestry_pcs_corrected.tsv (see step 3.5) -- they are NOT what
# gets passed to pcrelate() in step 4.
pcair_result <- pcair(
  gdsobj = genoData,
  kinobj = king_mat,
  divobj = king_mat
)

# Inspect how many samples went into the "unrelated" training set vs. the
# "related" set before trusting downstream results. pcrelate() (step 4)
# uses pcair_result$unrels as its training.set regardless of PC source.
summary(pcair_result)
cat(length(pcair_result$unrels), "samples in PC-AiR's unrelated set,",
    length(pcair_result$rels), "in the related set\n")

pcair_eigenvec <- as.data.frame(pcair_result$vectors)
colnames(pcair_eigenvec) <- paste0("PC", seq_len(ncol(pcair_eigenvec)))
pcair_eigenvec <- tibble(ID = rownames(pcair_result$vectors)) %>%
  bind_cols(pcair_eigenvec)
write_tsv(pcair_eigenvec, file.path(out_dir, "cohort_pcair.eigenvec"))

# Diagnostic-only scree info for PC-AiR's own PCs (see n_pcs_for_adjustment
# comment above for where to actually look to pick that value).
varprop_df <- tibble(
  PC = paste0("PC", seq_along(pcair_result$varprop)),
  varprop = pcair_result$varprop
)
write_tsv(varprop_df, file.path(out_dir, "cohort_pcair_varprop.txt"))
cat("PC-AiR variance proportion by PC (diagnostic only -- not used below):\n")
print(varprop_df)

# ---- Step 3.5: load the corrected, reference-projected ancestry PCs ----
# Written by r_scripts/ancestry_viz.R to ancestry/, alongside that step's
# other outputs. Sample IDs (IID) are matched and reordered against this
# GDS's own sample_ids -- not just assumed to line up -- and any mismatch
# fails loudly here rather than silently misaligning genotypes and PCs
# inside pcrelate().
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

# ---- Step 4: PC-Relate -- first pass, using step 3's automatic partition ----
genoIter <- GenotypeBlockIterator(genoData)

pcrelate_result_pass1 <- pcrelate(
  gdsobj = genoIter,
  pcs = ancestry_pcs_mat,
  training.set = pcair_result$unrels,
  BPPARAM = BiocParallel::SerialParam()
)
cat("pcrelate_result_pass1$kinBtwn columns:",
    paste(names(pcrelate_result_pass1$kinBtwn), collapse = ", "), "\n")

# ---- Step 4.5: convert pass-1 kinship back to matrix form for iteration ----
# PC-Relate's ancestry correction (via the `pcs` regression above) is
# applied to every pair, not just training-set pairs, so this is already a
# meaningfully better-corrected kinship estimate than raw KING even though
# pass 1's own training set (from step 3) is undersized/biased -- see the
# design-update comment near the top of this file.
kin_mat_pass1 <- pcrelateToMatrix(pcrelate_result_pass1, scaleKin = 2)

# ---- Step 5: diagnostic-only second pass, to flag confidently-related pairs ----
# Uses kin_mat_pass1 (ancestry-corrected) instead of raw king_mat for
# kinship, with GENESIS's own default kin.thresh/div.thresh (2^(-11/2)) --
# deliberately NOT loosened, unlike an earlier attempt at this step, which
# just admitted more true near-relatives into "unrelated" rather than
# fixing the actual ancestry-skew problem. This pass's own automatic
# partition may still be undersized for the same reason as pass 1 -- it's
# used here only to flag likely relatives for step 6, not as a final
# answer.
pcair_result_pass2 <- pcair(
  gdsobj = genoData,
  kinobj = kin_mat_pass1,
  divobj = king_mat
)
cat(length(pcair_result_pass2$unrels), "samples in pass 2's (still automatic) unrelated set,",
    length(pcair_result_pass2$rels), "in the related set\n")

genoIter2 <- GenotypeBlockIterator(genoData)
pcrelate_result_pass2 <- pcrelate(
  gdsobj = genoIter2,
  pcs = ancestry_pcs_mat,
  training.set = pcair_result_pass2$unrels,
  BPPARAM = BiocParallel::SerialParam()
)

kin_flagging <- pcrelate_result_pass2$kinBtwn %>%
  mutate(category = case_when(
    kin > 0.354  ~ "Duplicate/MZ twin",
    kin > 0.177  ~ "1st-degree",
    kin > 0.0884 ~ "2nd-degree",
    kin > 0.0442 ~ "3rd-degree",
    TRUE         ~ "Unrelated"
  ))
cat("Pass 2 (diagnostic) kinship categories, for flagging only:\n")
print(table(kin_flagging$category))

# ---- Step 6: manually build an ancestry-diverse unrelated training set ----
# related_thresh, n_training_clusters, n_per_cluster, and
# manual_force_include_ids are set near the top of this file.
flagged_ids <- kin_flagging %>%
  filter(kin > related_thresh) %>%
  select(ID1, ID2) %>%
  unlist() %>%
  unique()
cat(length(flagged_ids), "samples excluded from training-set candidacy",
    "(confidently related to someone, kin >", related_thresh, "):\n")
print(flagged_ids)

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
# just PC1/PC2) -- using km$centers directly rather than recomputing means
# also avoids any mismatch between the two.
candidates$dist_to_centroid <- sqrt(rowSums((scaled_mat - km$centers[km$cluster, ])^2))

training_selection <- candidates %>%
  group_by(cluster) %>%
  slice_min(order_by = dist_to_centroid, n = n_per_cluster) %>%
  ungroup()
cat(nrow(training_selection), "samples selected across", n_training_clusters, "clusters\n")
print(table(training_selection$cluster))

# Visual sanity check against the full candidate spread -- save rather
# than print(), since this script may run non-interactively via Rscript.
# Cross-check both plots against ancestry/ancestry_pc1_pc2.png and
# ancestry_pc3_pc4.png: if a visibly distinct ancestry group (e.g. an EAS
# singleton, AFR-leaning points) isn't covered by the red points below,
# add its sample ID(s) to manual_force_include_ids above and re-run from
# step 6.
p_train_pc12 <- ggplot() +
  geom_point(data = candidates, aes(PC1, PC2), color = "grey70", alpha = 0.5) +
  geom_point(data = training_selection, aes(PC1, PC2), color = "red", size = 3) +
  labs(title = "Manually-selected training set vs. candidate ancestry spread (PC1 vs PC2)") +
  theme_bw()
ggsave(file.path(out_dir, "training_set_pc1_pc2.png"), p_train_pc12, width = 8, height = 7, dpi = 150)

if (n_pcs_for_adjustment >= 4) {
  p_train_pc34 <- ggplot() +
    geom_point(data = candidates, aes(PC3, PC4), color = "grey70", alpha = 0.5) +
    geom_point(data = training_selection, aes(PC3, PC4), color = "red", size = 3) +
    labs(title = "Manually-selected training set vs. candidate ancestry spread (PC3 vs PC4)") +
    theme_bw()
  ggsave(file.path(out_dir, "training_set_pc3_pc4.png"), p_train_pc34, width = 8, height = 7, dpi = 150)
}

training_ids <- unique(c(training_selection$IID, manual_force_include_ids))
cat(length(training_ids), "final manually-selected training set:\n")
print(training_ids)
writeLines(training_ids, file.path(out_dir, "cohort_manual_training_set.txt"))

# ---- Step 7: final PC-AiR pass, forcing the manual set via unrel.set ----
# Confirmed against the GENESIS source (UW-GAC/GENESIS
# R/pcairPartition.R): unrel.set does NOT replace the automatic
# kin.thresh/div.thresh partition -- it forces the named samples into the
# unrelated set on top of it, and critically never re-flags two unrel.set
# members as "related to each other" even if their pairwise kinship (in
# kinobj) looks elevated, so an ancestry-inflated pair within the
# hand-picked set isn't second-guessed back out. GENESIS's own default
# thresholds are used here (no loosening needed anymore -- that's what
# unrel.set replaces).
pcair_result_final <- pcair(
  gdsobj = genoData,
  kinobj = kin_mat_pass1,
  divobj = king_mat,
  unrel.set = training_ids
)
summary(pcair_result_final)
cat(length(pcair_result_final$unrels), "samples in the FINAL unrelated training set,",
    length(pcair_result_final$rels), "in the related set\n")

# ---- Step 8: final PC-Relate pass -- this is the result actually written out ----
genoIter3 <- GenotypeBlockIterator(genoData)

pcrelate_result_final <- pcrelate(
  gdsobj = genoIter3,
  pcs = ancestry_pcs_mat,
  training.set = pcair_result_final$unrels,
  BPPARAM = BiocParallel::SerialParam()
)

kin_adjusted <- pcrelate_result_final$kinBtwn %>%
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
