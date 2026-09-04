# Run interactively (e.g. `module load R/4.4.0 && R`), NOT via sbatch or
# Rscript in a batch job. Bioconductor installs can take a while and
# sometimes prompt to update dependent packages ("Update all/some/none?
# [a/s/n]:") -- that prompt hangs forever in a non-interactive job.
# One-time setup for jobs/genesis_pcrelate.sh; packages install into your
# personal R library.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("GENESIS", "GWASTools", "SNPRelate", "gdsfmt"))

# dplyr/readr (tidyverse) are assumed already available, since
# r_scripts/relatedness_viz.R and crosscheckfingerprint_scores_viz.R
# already depend on tidyverse.
