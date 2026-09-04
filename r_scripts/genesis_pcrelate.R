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

# ---- Step 4: PC-Relate -- ancestry-adjusted kinship ----
genoIter <- GenotypeBlockIterator(genoData)

pcrelate_result_1 <- pcrelate(
  gdsobj = genoIter,
  pcs = pcair_result$vector[, 1:n_pcs_for_adjustment],
  training.set = pcair_result$unrels,
  BPPARAM = BiocParallel::SerialParam()
)

# ---- Step 4.5: iterate -- use this first-pass kinship to refine PC-AiR ----
kin_mat_1 <- pcrelateToMatrix(pcrelate_result_1, scaleKin = 2)  # convert kinBtwn output back to matrix form

pcair_result_3 <- pcair(
  gdsobj = genoData,
  kinobj = kin_mat_1,
  divobj = king_mat,
  kin.thresh = 2^(-7/2),   # ~0.0442, conventional 3rd-degree cutoff, instead of the 4th-degree default
  div.thresh = -2^(-7/2)
)

summary(pcair_result_3)  # check unrelated-set size now

# ---- Step 5: second, refined PC-Relate pass ----
genoIter2 <- GenotypeBlockIterator(genoData)

pcrelate_result_3 <- pcrelate(
  gdsobj = genoIter2,
  pcs = pcair_result_3$vectors[, 1:n_pcs_for_adjustment],
  training.set = pcair_result_3$unrels,
  BPPARAM = BiocParallel::SerialParam()
)

kin_adjusted <- pcrelate_result_3$kinBtwn %>%
  mutate(category = case_when(
    kin > 0.354  ~ "Duplicate/MZ twin",
    kin > 0.177  ~ "1st-degree",
    kin > 0.0884 ~ "2nd-degree",
    kin > 0.0442 ~ "3rd-degree",
    TRUE         ~ "Unrelated"
  ))

table(kin_adjusted$category)

ggplot(kin_adjusted, aes(k0, kin)) +
  geom_hline(yintercept = 2^(-seq(3,9,2)/2), 
             linetype = "dashed", 
             color = "grey") +
  geom_point(alpha = 0.5) +
  theme_bw()

print(table(kin_adjusted$category))

write_tsv(kin_adjusted, out_fn)
cat("Wrote", nrow(kin_adjusted), "pairwise estimates to", out_fn, "\n")

close(gds_reader)
