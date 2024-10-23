#################################
# Author: Xiangjie Zhao
# Date: 2023.04.06
#################################
rm(list = ls())
set.seed(1)
library(tidyverse)
setwd('~/project/hormone_receptor_scRNA-2022-12-11/')

#################################
# FDA-approved drug targets
df_FDA <- read_tsv('input/protein_class_FDA.tsv', col_names = T)
vec_FDA <- df_FDA$Gene
length(vec_FDA)

#################################
# HR list and Non-HR list
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
vec_receptor <- df_receptor$Gene_symbol
length(vec_receptor)
df_receptor$Gene_info <- 'HR'

df_all <- readxl::read_excel('input/reference-genome/proteinCodingGeneSymbol.xlsx', col_names = F)
colnames(df_all) <- 'Gene_symbol'
df <- left_join(df_all, df_receptor,"Gene_symbol")
df <- df %>% replace_na(list(Gene_info = 'Non-HR'))
df$`Gene_info` <- factor(df$`Gene_info`, levels = c('Non-HR','HR'))
df <- df %>% mutate(`FDA-approved targets` = if_else(Gene_symbol %in% vec_FDA, "Yes", "No"))
df$`FDA-approved targets` <- factor(df$`FDA-approved targets`, levels = c('No','Yes'))
df$`Receptor_subcellular` <- factor(df$`Receptor_subcellular`, levels = c('Nucleus','Membrane'))

ggplot(df, aes(Gene_info, fill=`FDA-approved targets`)) +
  geom_bar(position = 'fill', width = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  scale_fill_manual(values = c('#999999','#59BBCC')) +
  coord_flip() +
  labs(x='Protein-coding genes', y='Percentage')
ggsave('output/HPA-single_cell/figure1/FDA-allGene.pdf', width = 4.5, height = 1.5)

#################################
# Memebrane HR and Nucleus HR
df <- df %>% filter(Gene_info == 'HR')
ggplot(df, aes(Receptor_subcellular, fill=`FDA-approved targets`)) +
  geom_bar(position = 'fill', width = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  scale_fill_manual(values = c('#999999','#59BBCC')) +
  coord_flip() +
  labs(x='HR types', y='Percentage')
ggsave('output/HPA-single_cell/figure1/FDA-HRtypes.pdf', width = 4.5, height = 1.5)

