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
# NOT recognize plink2's --make-king-table column names (IID1, IID2,
# KINSHIP) -- despite what some tutorials assume, it does not autodetect
# or support the plink2 format directly (confirmed against the GENESIS
# source, UW-GAC/GENESIS R/makeSparseMatrix.R: it does a strict
# intersect() against literal "ID1"/"ID2"/<estimator> column names). This
# script renames the columns to what kingToMatrix expects before calling
# it, rather than assuming compatibility.

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
# GENESIS vignette's own example uses 2; adjust this after looking at
# ref_pca/cohort_projected_pca (or the PC-AiR PCs computed below, in
# genesis/cohort_pcair.eigenvec) to see how many PCs actually separate
# distinct ancestry clusters in this cohort -- there's no way to pick
# this correctly without a human looking at the plot.
n_pcs_for_adjustment <- 2

# ---- Step 1: convert the pruned, QC'd cohort genotypes to GDS format ----
# bed_prefix.bed/.bim/.fam is written by jobs/genesis_pcrelate.sh via
# `plink2 --pfile relatedness/cohort_qc --extract
# relatedness/cohort_pruned.prune.in --make-bed` -- the same pruned
# marker set already used to compute cohort_king.kin0, so the preliminary
# KING estimate and this re-analysis are on consistent footing.
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
# Rename plink2's --make-king-table columns to what kingToMatrix expects.
# matches()-based renaming (regex, not exact-string) so this is robust to
# plink2's leading "#" on the first header column, without needing to
# confirm that convention ahead of time. read_table() (whitespace-
# flexible), not read_tsv(), matching how r_scripts/relatedness_viz.R
# already successfully reads this exact file.
king_raw <- read_table("relatedness/cohort_king.kin0", show_col_types = FALSE)
cat("cohort_king.kin0 columns:", paste(names(king_raw), collapse = ", "), "\n")

king_renamed <- king_raw %>%
  rename(
    ID1 = matches("IID1$"),
    ID2 = matches("IID2$"),
    Kinship = KINSHIP
  )
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
