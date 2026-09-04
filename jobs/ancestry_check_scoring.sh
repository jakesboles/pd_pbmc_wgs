#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name ancestry_check_scoring
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 32G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A.log
#SBATCH --verbose

set -euo pipefail

cd /projects/b1169/boles/pd_pbmc_wgs

module load plink/2.001

plink2 \
  --threads 8 \
  --pfile ancestry/ref_shared \
  --extract ancestry/ref_pruned.prune.in \
  --read-freq ancestry/ref_pca.afreq \
  --score ancestry/ref_pca.eigenvec.var 2 3 header-read no-mean-imputation variance-standardize \
  --score-col-nums 5-14 \
  --out ancestry/ref_selfprojected_pca