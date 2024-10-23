#################################
# Author: Xiangjie Zhao
# Date: 2023.03.24
# Title: Analysis of hormone receptor gene list
#################################

rm(list = ls())
set.seed(1)
library(tidyverse)
library(RColorBrewer)
library(gganatogram)
library(dplyr)
library(viridis)
library(gridExtra)
setwd('~/project/hormone_receptor_scRNA-2022-12-11/')

#################################
gtex_hr <- read_tsv('input/gtex/GTEx_HR_combined.tsv')
gtex_hr <- gtex_hr %>% group_by(Description,SMTS,SEX) %>% filter(TPM > (median(TPM) - sd(TPM)), TPM < median(TPM) + sd(TPM)) # 去离群值
gtex_hr$organ <- gtex_hr$SMTS
gtex_hr$organ <- gtex_hr$organ %>% recode(`Adipose Tissue` = 'adipose_tissue',
                                          `Adrenal Gland` = 'adrenal_gland',
                                          `Bladder` = 'urinary_bladder',
                                          `Blood` = 'leukocyte',
                                          `Blood Vessel` = 'aorta',
                                          `Brain` = 'brain',
                                          `Breast` = 'breast',
                                          `Colon` = 'colon',
                                          `Esophagus` = 'esophagus',
                                          `Heart` = 'heart',
                                          `Kidney` = 'kidney',
                                          `Liver` = 'liver',
                                          `Lung` = 'lung',
                                          `Muscle` = 'skeletal_muscle',
                                          `Nerve` = 'nerve',
                                          `Pancreas` = 'pancreas',
                                          `Pituitary` = 'pituitary_gland',
                                          `Prostate` = 'prostate',
                                          `Salivary Gland` = 'salivary_gland',
                                          `Skin` = 'skin',
                                          `Small Intestine` = 'small_intestine',
                                          `Spleen` = 'spleen',
                                          `Stomach` = 'stomach',
                                          `Testis` = 'testis',
                                          `Thyroid` = 'thyroid_gland')

#################################
# compute mean expression in tissues
gene <- 'NR1H2'
sex <- 'Male'
gtex_hr_gene <- gtex_hr %>% filter(Description == gene, SEX == sex) %>% group_by(organ) %>% summarize(TPM = mean(TPM))
gtex_hr_gene <- gtex_hr_gene %>% arrange(TPM)
gtex_hr_gene <- rename(gtex_hr_gene, value = TPM)
# gtex_hr_gene$value <- log10(gtex_hr_gene$value)

gganatogram(data=gtex_hr_gene, fillOutline='#a6bddb', organism='human', sex='male', fill="value") + 
  theme_void() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title=gene) +
  theme(legend.position='none', plot.title = element_text(hjust = 0.5, vjust=-7))

ggsave(paste0('output/HPA-single_cell/figure3/gganatogram-',sex,'-',gene,'.pdf'), width = 5, height = 8)
# 
# 
# ################################
# # tissue color
# # 74 colors
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# gtex_hr_gene$colour <- col_vector[1:dim(gtex_hr_gene)[1]]
# gtex_hr_gene <- gtex_hr_gene %>% arrange(organ)
# gtex_hr_gene$organ <- factor(gtex_hr_gene$organ, levels=gtex_hr_gene$organ)
# gganatogram(data=gtex_hr_gene, fillOutline='white', organism='human', sex='male', fill='colour') + theme_void()
# ggsave('output/HPA-single_cell/figure3/gganatogram-tissue.pdf', width = 5, height = 7)
# 
# ggplot(gtex_hr_gene, aes(organ,organ, color=organ))+
#   geom_point()+
#   scale_color_manual(values = gtex_hr_gene$colour) +
#   theme_bw()
# ggsave('output/HPA-single_cell/figure3/gganatogram-tissue-legend.pdf', width = 5, height = 8)
