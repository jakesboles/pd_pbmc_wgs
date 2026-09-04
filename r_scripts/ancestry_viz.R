library(tidyverse)
library(ggplot2)

setwd("/projects/b1169/boles/pd_pbmc_wgs")

# compute correction factor from plink2's scoring

ref_eigenvec <- read_tsv("ancestry/ref_pca.eigenvec")          # ground truth: PLINK2's own --pca output
ref_selfproj <- read_tsv("ancestry/ref_selfprojected_pca.sscore")  # same individuals, scored via --score
eigenval     <- read_tsv("ancestry/ref_pca.eigenval", col_names = "eigenval")

merged <- ref_selfproj %>%
  rename(IID = `#IID`) %>%
  inner_join(
    ref_eigenvec %>% rename(IID = `#IID`),
    by = "IID",
    suffix = c("_score", "_eigenvec")
  )

n_pcs <- 10
correction <- tibble(
  PC = paste0("PC", seq_len(n_pcs)),
  slope = NA_real_,
  r_squared = NA_real_,
  theoretical = -2 / sqrt(eigenval$eigenval[seq_len(n_pcs)])
)

for (k in seq_len(n_pcs)) {
  score_col <- paste0("PC", k, "_AVG")
  eig_col   <- paste0("PC", k)
  
  fit <- lm(merged[[eig_col]] ~ merged[[score_col]])
  
  correction$slope[k]     <- coef(fit)[2]       # the slope IS the correction factor
  correction$r_squared[k] <- summary(fit)$r.squared
}

print(correction)

# ---- Confirm the fit is essentially perfect before trusting the correction ----
if (any(correction$r_squared < 0.999)) {
  warning("At least one PC's self-projection didn't recover its own eigenvector cleanly — inspect before proceeding.")
}

correction_factor <- correction$slope

# apply correction factor and visualize

# ---- Load inputs ----
cohort   <- read_tsv("ancestry/cohort_projected_pca.sscore")

apply_correction <- function(df, n_pcs = 10) {
  for (k in seq_len(n_pcs)) {
    avg_col <- paste0("PC", k, "_AVG")
    df[[paste0("PC", k)]] <- df[[avg_col]] * correction_factor[k]
  }
  df
}

ref_corrected    <- apply_correction(ref_selfproj)
cohort_corrected <- apply_correction(cohort)

# ---- Combine into one long data frame for plotting ----
ref_plot <- ref_corrected %>%
  dplyr::rename(IID = `#IID`) %>%
  select(IID, SuperPop, PC1, PC2, PC3, PC4) %>%
  mutate(group = SuperPop)

cohort_plot <- cohort_corrected %>%
  # dplyr::rename(IID = `#IID`) %>%
  select(IID, PC1, PC2, PC3, PC4) %>%
  mutate(group = "Cohort")

combined <- bind_rows(ref_plot, cohort_plot)

# ---- PC1 vs PC2 ----
p1 <- ggplot() +
  geom_point(data = filter(combined, group != "Cohort"),
             aes(x = PC1, y = PC2, color = group),
             alpha = 0.6, size = 2) +
  geom_point(data = filter(combined, group == "Cohort"),
             aes(x = PC1, y = PC2),
             color = "black", shape = 4, size = 2, stroke = 1) +
  labs(color = "1000G SuperPop") +
  ggtitle("PC1 vs PC2") +
  theme_bw() + 
  theme(axis.text = element_text(color = "black"))

ggsave("ancestry/ancestry_pc1_pc2.png", p1, width = 8, height = 7, dpi = 150)

# ---- PC3 vs PC4, same treatment ----
p2 <- ggplot() +
  geom_point(data = filter(combined, group != "Cohort"),
             aes(x = PC3, y = PC4, color = group),
             alpha = 0.6, size = 2) +
  geom_point(data = filter(combined, group == "Cohort"),
             aes(x = PC3, y = PC4),
             color = "black", shape = 4, size = 2, stroke = 1) +
  labs(color = "1000G SuperPop", title = "PC3 vs PC4") +
  ggtitle("PC3 vs PC4") +
  theme_bw() + 
  theme(axis.text = element_text(color = "black"))

ggsave("ancestry/ancestry_pc3_pc4.png", p2, width = 8, height = 7, dpi = 150)

# ---- Save the corrected, labeled cohort table for downstream use (e.g. QTL covariates) ----
write_tsv(cohort_corrected %>% select(IID, PC1:PC10),
          "cohort_ancestry_pcs_corrected.tsv")

print(p1)