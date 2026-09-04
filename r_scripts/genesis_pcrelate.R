# Ancestry-adjusted relatedness via GENESIS's PC-AiR / PC-Relate pipeline.
# Run by jobs/genesis_pcrelate.sh. Requires the packages installed by
# r_scripts/install_genesis_packages.R (run that first, interactively).
#
# Note on "ancestry" here: PC-AiR computes its OWN ancestry-representative
# PCs directly from the pruned cohort genotypes below -- it does NOT reuse
# ancestry/cohort_projected_pca.sscore (the 1000G-projected PCs from
# ancestry_pca.sh). That's standard/correct GENESIS methodology, not an
# oversight: PC-AiR needs PCs computed consistently with its own kinship
# matrix to identify a maximally-unrelated "training" subset, which a
# separately-computed, differently-pruned set of projected PCs wouldn't
# guarantee. This step is "ancestry-aware" relatedness estimation in the
# sense GENESIS defines it, not a literal reuse of the earlier PCA output.
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

# How many PC-AiR PCs to hand to PC-Relate for ancestry adjustment. The
# GENESIS vignette's own example uses 2, but that's not a default to trust
# blindly -- pick this by looking at two things, both written out below
# before pcrelate() runs:
#   1. genesis/cohort_pcair_varprop.txt -- each PC's proportion of variance
#      explained (pcair_result$varprop). Look for the "elbow" where added
#      PCs stop explaining much more variance -- PCs past that point are
#      mostly noise, not structure.
#   2. genesis/cohort_pcair.eigenvec -- pairs-plot PC1 vs PC2, PC3 vs PC4,
#      etc. (e.g. with r_scripts/relatedness_viz.R-style ggplot code) and
#      look for how many PCs still visibly separate distinct clusters,
#      the same way you already did for ancestry_pca.sh's 1000G-projected
#      PCs against known SuperPop labels.
# Confirmed against the GENESIS source (UW-GAC/GENESIS R/pcair.R): the
# object pcair() returns has $values (eigenvalues) and $varprop (proportion
# of variance explained) fields, in addition to $vectors/$unrels/$rels
# used elsewhere in this script.
n_pcs_for_adjustment <- 2

# ---- Step 1: convert the pruned, QC'd cohort genotypes to GDS format ----
# bed_prefix.bed/.bim/.fam is written by jobs/genesis_pcrelate.sh via
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

king_mat <- kingToMatrix(
  king_renamed_fn,
  estimator = "Kinship",
  sample.include = sample_ids
)

# ---- Step 3: PC-AiR -- ancestry PCs robust to relatedness ----
# Uses the same KING matrix for both kinship AND divergence, per GENESIS
# convention: KING-robust kinship already encodes ancestry divergence in
# its negative values.
pcair_result <- pcair(
  gdsobj = genoData,
  kinobj = king_mat,
  divobj = king_mat
)

# Inspect how many samples went into the "unrelated" training set vs. the
# "related" set before trusting downstream results.
summary(pcair_result)
cat(length(pcair_result$unrels), "samples in PC-AiR's unrelated set,",
    length(pcair_result$rels), "in the related set\n")

pcair_eigenvec <- as.data.frame(pcair_result$vectors)
colnames(pcair_eigenvec) <- paste0("PC", seq_len(ncol(pcair_eigenvec)))
pcair_eigenvec <- tibble(ID = rownames(pcair_result$vectors)) %>%
  bind_cols(pcair_eigenvec)
write_tsv(pcair_eigenvec, file.path(out_dir, "cohort_pcair.eigenvec"))

# Scree diagnostic for picking n_pcs_for_adjustment above -- see that
# variable's comment. Written before pcrelate() runs so it's available to
# review even if pcrelate() itself fails.
varprop_df <- tibble(
  PC = paste0("PC", seq_along(pcair_result$varprop)),
  varprop = pcair_result$varprop
)
write_tsv(varprop_df, file.path(out_dir, "cohort_pcair_varprop.txt"))
cat("PC-AiR variance proportion by PC:\n")
print(varprop_df)

# ---- Step 4: PC-Relate -- ancestry-adjusted kinship ----
genoIter <- GenotypeBlockIterator(genoData)

pcrelate_result <- pcrelate(
  gdsobj = genoIter,
  pcs = pcair_result$vectors[, seq_len(n_pcs_for_adjustment), drop = FALSE],
  training.set = pcair_result$unrels,
  BPPARAM = BiocParallel::SerialParam()
)

# ---- Step 5: extract and categorize pairwise kinship, same thresholds as plink_relatedness.sh ----
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

write_tsv(kin_adjusted, out_fn)
cat("Wrote", nrow(kin_adjusted), "pairwise estimates to", out_fn, "\n")

close(gds_reader)
