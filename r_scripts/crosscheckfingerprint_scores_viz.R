library(tidyverse)
library(ggplot2)

setwd("/projects/b1169/boles/pd_pbmc_wgs")

files <- list.files("crosscheck",
                    recursive = F, full.names = T,
                    pattern = "crosscheck_metrics")

df <- map(files,
          read_table,
          skip = 6)

df <- df %>% 
  list_rbind()

write.csv(df,
          file = "crosscheck/crosscheck_metrics_compiled.csv",
          row.names = F)

hist(df$LOD_SCORE)

table(df$RESULT)

df %>% 
  ggplot(aes(x = LOD_SCORE)) + 
  geom_histogram(color = "black", fill = "cadetblue",
                 binwidth = 2) + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(15, 90, 5)) +
  labs(y = "N",
       x = "LOD score") +
  theme_linedraw()
ggsave(filename = "crosscheck/lod_histogram.png",
       units = "in", dpi = 600,
       height = 4, width = 6)
  