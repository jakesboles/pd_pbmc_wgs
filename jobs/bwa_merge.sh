#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name bwa_merge
#SBATCH --nodes 1
#SBATCH --array=59,99
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/pd_pbmc_wgs/logs/%x_%A_%a.log
#SBATCH --verbose

cd /projects/b1169/boles/pd_pbmc_wgs

module load samtools/1.16.1-gcc-10.4.0
module load gatk/4.4.0.0

PARAMS_FILE="bowtie_params_id.txt"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f1 -d,)

files=$(ls "$PWD"/bwa_bam/* | grep "${sample}[_]" | grep ".bai" -Ev)

echo ${sample}
echo ${files}

echo "Merging BAMs"

samtools merge \
-@ 16 \
bwa_bam/${sample}.merged.bam \
${files} \

echo "Indexing merged BAM"

samtools index bwa_bam/${sample}.merged.bam

echo "Sorting merged BAM on coordinate"

samtools sort -@ 16 \
  -o bwa_bam/${sample}.merged.sorted.bam \
  bwa_bam/${sample}.merged.bam
  
echo "Completed successfully"

