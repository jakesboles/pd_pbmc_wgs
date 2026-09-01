#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name genomicsdbiimport
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --array=1-23
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

module load gatk/4.4.0.0

cd /projects/b1169/boles/pd_pbmc_wgs

CHROMOSOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "chromosomes.txt")

gatk --java-options "-Xmx32g -Xms32g" GenomicsDBImport \
  --sample-name-map cohort.sample_map \
  --genomicsdb-workspace-path "genomics_db/chr${CHROMOSOME}_db" \
  -L "chr${CHROMOSOME}" \
  --reader-threads 16