cd /projects/b1042/Gate_Lab/boles/pd_wgs

ls haplotype_caller/*.g.vcf.gz | while read gvcf; do   
  sample=$(basename "$gvcf" .output.g.vcf.gz);   
  echo -e "${sample}\t${gvcf}"; 
done > cohort.sample_map