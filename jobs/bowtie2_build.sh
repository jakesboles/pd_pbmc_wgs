#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name bowtie2_build
#SBATCH --nodes 1
#SBATCH --array=1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1042/Gate_Lab/boles/pd_wgs/logs/%x_%A.log
#SBATCH --verbose

cd /projects/p31535/boles/bowtie2_ref_v2

module load bowtie2/2.5.4

bowtie2-build \
  --threads 16 \
  ../Homo_sapiens_assembly38.fasta \
  Homo_sapiens_assembly38
