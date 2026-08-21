library(tidyverse)

setwd("/projects/b1042/Gate_Lab/boles/pd_wgs")

files <- list.files("fastq")

# bowtie2 params -----------------------------------------------------------

ids <- str_split_i(files, "_", i = 1) %>% unique()

mat <- matrix(nrow = length(ids),
              ncol = 3)

for (i in seq_along(ids)){
  all_reads <- files[str_detect(files, ids[i])]
  
  mat[i, 1] <- ids[i]
  
  mat[i, 2] <- paste(all_reads[str_detect(all_reads, "_R1_")],
                     collapse = ",")
  mat[i, 3] <- paste(all_reads[str_detect(all_reads, "_R2_")],
                     collapse = ",")
}

mat[, 1] %>% 
  as_tibble() %>% 
  write.table("bowtie_params_id.txt",
              quote = F,
              row.names = F,
              col.names = F,
              sep = "|")
mat[, 2] %>% 
  as_tibble() %>% 
  write.table("bowtie_params_r1.txt",
              quote = F,
              row.names = F,
              col.names = F,
              sep = "|")
mat[, 3] %>% 
  as_tibble() %>% 
  write.table("bowtie_params_r2.txt",
              quote = F,
              row.names = F,
              col.names = F,
              sep = "|")

# cutadapt params ---------------------------------------------------------

reps <- str_remove_all(files, "_R1_001.fastq.gz") %>%
  str_remove_all("_R2_001.fastq.gz") %>%
  unique()

mat <- matrix(nrow = length(reps),
              ncol = 3)

for (i in seq_along(reps)) {
  fastq_pair <- files[str_detect(files, reps[i])]
  
  mat[i, 1] <- fastq_pair[str_detect(fastq_pair, "R1")]
  mat[i, 2] <- fastq_pair[str_detect(fastq_pair, "R2")]
  mat[i, 3] <- reps[i]

}

mat %>%
  as_tibble() %>% 
  write.table(file = "cutadapt_params.txt",
              quote = F,
              row.names = F,
              col.names = F,
              sep = ",")

# BWA params --------------------------------------------------------------

as.data.frame(mat) %>% 
  mutate(lane = str_split_i(V3, "_", i = 3),
         sample = str_split_i(V3, "_", i = 1)) %>%
  write.table(file = "bwa_params.txt",
              quote = F,
              row.names = F,
              col.names = F,
              sep = ",")
