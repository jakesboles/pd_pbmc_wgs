cd /projects/b1169/boles/pd_pbmc_wgs

mkdir -p params

ls haplotype_caller/*.g.vcf.gz | while read gvcf; do
  sample=$(basename "$gvcf" .output.g.vcf.gz);
  echo -e "${sample}\t${gvcf}";
done > params/cohort.sample_map