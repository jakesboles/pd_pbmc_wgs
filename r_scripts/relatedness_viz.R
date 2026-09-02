library(tidyverse)
library(ggplot2)

setwd("/projects/b1169/boles/pd_pbmc_wgs")

df <- read_table("relatedness/cohort_king.kin0")

head(df)

ggplot(df,
       aes(x = KINSHIP,
           y = IBS0)) + 
  geom_point() + 
  theme_linedraw()
