#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name bwa
#SBATCH --nodes 1
#SBATCH --array=1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1042/Gate_Lab/boles/pd_wgs/logs/%x_%A.log
#SBATCH --verbose

cd /projects/p31535/boles

module load bwa/0.7.17
module load samtools/1.16.1-gcc-10.4.0
module load gatk/4.4.0.0

bwa-mem2-2.2.1_x64-linux/bwa-mem2 index Homo_sapiens_assembly38.fasta

samtools faidx Homo_sapiens_assembly38.fasta

gatk CreateSequenceDictionary -R Homo_sapiens_assembly38.fasta
