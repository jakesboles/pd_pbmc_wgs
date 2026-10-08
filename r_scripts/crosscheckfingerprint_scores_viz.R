library(tidyverse)
library(ggplot2)

setwd("/projects/b1169/boles/pd_pbmc_wgs")

# df <- read_table("crosscheck/cohort_all_pairs.crosscheck_metrics",
#                  skip = 6)
df <- read_table("crosscheck/cohort_all_pairs_gex.crosscheck_metrics",
                 skip = 6)

# write.csv(df,
#           file = "crosscheck/crosscheck_metrics_compiled.csv",
#           row.names = F)

hist(df$LOD_SCORE)

table(df$RESULT)

df %>% 
  ggplot(aes(x = LEFT_GROUP_VALUE,
             y = RIGHT_GROUP_VALUE)) +
  geom_tile(aes(fill = LOD_SCORE)) + 
  scale_fill_gradient2(high = "midnightblue") +
  labs(y = "ATAC library", # y = "GEX library"
       x = "WGS library") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
# ggsave(filename = "crosscheck/lod_heatmap.png",
#        units = "in", dpi = 600,
#        height = 4, width = 6)
ggsave(filename = "crosscheck/lod_gex_heatmap.png",
       units = "in", dpi = 600,
       height = 4, width = 6)

df %>% 
  ggplot(aes(x = LOD_SCORE)) + 
  geom_histogram(color = "black", fill = "cadetblue",
                 binwidth = 2) + 
  scale_y_continuous(expand = c(0, 0)) +
  # scale_x_continuous(breaks = seq(15, 90, 5)) +
  labs(y = "N",
       x = "LOD score") +
  theme_linedraw()
# ggsave(filename = "crosscheck/lod_histogram.png",
#        units = "in", dpi = 600,
#        height = 4, width = 6)
ggsave(filename = "crosscheck/lod_gex_histogram.png",
       units = "in", dpi = 600,
       height = 4, width = 6)
  