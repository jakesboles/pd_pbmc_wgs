#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name genotype_gvcfs
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --array=1-23
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

module load gatk/4.4.0.0

cd /projects/b1169/boles/pd_pbmc_wgs

CHROMOSOME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "params/chromosomes.txt")
DB="gendb://genomics_db/chr${CHROMOSOME}_db"

 gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /projects/p31535/boles/Homo_sapiens_assembly38.fasta \
   -V ${DB} \
   -O genotyped_gvcfs/chr${CHROMOSOME}.vcf.gz