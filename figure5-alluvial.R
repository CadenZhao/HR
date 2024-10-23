#################################
# Author: Xiangjie Zhao
# Date: 2023.03.24
# Title: Analysis of hormone receptor gene list
#################################

rm(list = ls())
set.seed(1)
library(tidyverse)
library(ggalluvial)
setwd('~/project/hormone_receptor_scRNA-2022-12-11/')

#################################
# 读取激素和HR的信息
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
df_receptor$Hormone_name1 <- recode(df_receptor$Hormone_name1, 'orphan'='Unknown')
df_receptor$Hormone_chemical_classes <- recode(df_receptor$Hormone_chemical_classes, 'NA'='Unknown')
df_receptor$Hormone1_tissue <- recode(df_receptor$Hormone1_tissue, 'NA'='Unknown')

# alluvial plot 
df_receptor <- df_receptor %>% separate_rows(Hormone1_tissue, sep = ',')
df_receptor <- df_receptor %>% arrange(Hormone_chemical_classes,Hormone1_tissue)
df_receptor$Hormone_chemical_classes <- factor(df_receptor$Hormone_chemical_classes, levels = unique(df_receptor$Hormone_chemical_classes))
ggplot(data = df_receptor, aes(axis1 = Hormone_name1, axis2 = Hormone_chemical_classes, axis3=Hormone1_tissue, axis4 = Receptor_subcellular, axis5 = Gene_symbol)) +
  geom_alluvium(aes(fill = Receptor_subcellular)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void()

ggsave('output/HPA-single_cell/figure5/alluvial.pdf', width = 15, height = 15)
