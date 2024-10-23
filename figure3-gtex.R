#################################
# Author: Xiangjie Zhao
# Date: 2023.03.25
# Title: Analysis of hormone receptor gene list
#################################

rm(list = ls())
set.seed(1)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggpubr)
library(comprehenr)
setwd('~/project/hormone_receptor_scRNA-2022-12-11/')
color=rev(brewer.pal(11, name = 'RdYlBu')) # define a heatmap color

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

common_tissue <- c("Thyroid","Stomach","Spleen","Small Intestine","Skin","Salivary Gland","Pituitary","Pancreas","Nerve","Muscle","Lung","Liver","Kidney","Heart","Esophagus","Colon","Breast","Brain","Blood Vessel","Blood","Bladder","Adrenal Gland","Adipose Tissue")

#################################
gtex_hr <- read_tsv('input/gtex/GTEx_HR_combined.tsv')
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
df_receptor$Hormone_chemical_classes <- recode(df_receptor$Hormone_chemical_classes, 'NA'='Orphan')
df_receptor$Description <- df_receptor$Gene_symbol
gtex_hr <- left_join(gtex_hr, df_receptor, by='Description')

################################
# 设置阈值，基因在组织中的平均TPM低于此值视为不表达
df <- gtex_hr %>% group_by(Description, SMTS) %>% summarise(TPM=mean(TPM))
# 这里即log1p(TPM) < 0.1视为不表达
ggplot(df, aes(log1p(TPM)))+
  geom_density()+
  geom_vline(xintercept=0.1)+
  labs(y='Distribution density')+
  theme_classic()
# ggsave('output/HPA-single_cell/figure3/gtex-cutoff.pdf', width = 2.5, height = 2)

#################################
# 30种组织表达HR数目的排名（分HR核膜类型）
df <- gtex_hr %>% group_by(Description, SMTS, Receptor_subcellular) %>% summarise(TPM=mean(TPM)) # 计算HR在组织（SMTS）中的表达平均值
df <- df %>% group_by(SMTS, Receptor_subcellular) %>% summarise(count = sum(log1p(TPM) > 0.1)) # 计算组织所表达的HR数目
level <- unique(arrange(df %>% group_by(SMTS) %>% summarise(count = sum(count)), count)$SMTS)
df$SMTS <- factor(df$SMTS, levels = level)
ggplot(df, aes(SMTS, count, fill=Receptor_subcellular)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  # coord_flip() +
  labs(x='Tissue', y='No. of HRs', fill='HR types') +
  scale_y_continuous(breaks = seq(0,140,10)) +
  scale_fill_manual(values = c("#ED1B28","#237EB6")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure3/gtex-number-of-expressed-receptors-HRtype-1.pdf', width = 7, height = 5)

# 132 HRs expressing in testis
HR_testis <- df %>% filter(SMTS == 'Testis', log1p(TPM) > 0.1) %>% pull(Description)
setdiff(df_receptor$Gene_symbol, HR_testis) # 6 HRs not expressing in testis: [1] "AGTR2"  "AVPR1B" "GHSR"   "NPY4R"  "RXFP3"  "TRHR"  

#################################
# 30种组织表达HR数目的排名（分性别）
df <- gtex_hr %>% group_by(Description, SMTS, SEX) %>% summarise(TPM=mean(TPM)) # 计算HR在组织（SMTS）中的表达平均值
df <- df %>% group_by(SMTS, SEX) %>% summarise(count = sum(log1p(TPM) > 0.1)) # 计算组织所表达的HR数目
level <- unique(arrange(df %>% group_by(SMTS) %>% summarise(count = sum(count)), count)$SMTS)
df$SMTS <- factor(df$SMTS, levels = level)
ggplot(df, aes(SMTS, count, fill=SEX)) +
  geom_bar(stat = 'identity', position = "dodge", width = 0.7) +
  theme_classic() +
  labs(x='Tissue', y='Number of expressed HRs', fill='SEX') +
  scale_fill_manual(values = c("#3BAF51","#9C4FA1")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave('output/HPA-single_cell/figure3/gtex-number-of-expressed-receptors-SEX.pdf', width = 8, height = 4)

df_female <- df %>% filter(SEX == 'Female', SMTS %in% common_tissue)
df_male <- df %>% filter(SEX == 'Male', SMTS %in% common_tissue)
df_female$change_proportion <- (df_male$count - df_female$count) / df_female$count * 100 # 男性相比女性各种组织所表达HR数目变化的百分比（负值代表减少）
df_female$change_direction <- to_vec(for(i in df_female$change_proportion) if(i<0) 'Down' else 'Up')
ggplot(df_female, aes(SMTS, change_proportion, fill=change_direction)) +
  geom_bar(stat = 'identity', color='black', linewidth=0.3) +
  coord_flip() +
  theme_classic() +
  labs(x='Tissue', y='Changed percentage', fill='Change direction') +
  scale_fill_manual(values = c('#237EB5','#EE1B27'))
# ggsave('output/HPA-single_cell/figure3/gtex-number-of-expressed-receptors-SEX-changePercentageMale2Female.pdf', width = 4, height = 4)

#################################
# 30种组织表达HR数目的排名（分年龄）
df <- gtex_hr %>% group_by(Description, SMTS, AGE) %>% summarise(TPM=mean(TPM)) # 计算HR在组织（SMTS）中的表达平均值
df <- df %>% group_by(SMTS, AGE) %>% summarise(count = sum(log1p(TPM) > 0.1)) # 计算组织所表达的HR数目
level <- unique(arrange(df %>% group_by(SMTS) %>% summarise(count = sum(count)), count)$SMTS)
df$SMTS <- factor(df$SMTS, levels = level)
ggplot(df, aes(SMTS, count, fill=AGE)) +
  geom_bar(stat = 'identity', position = "dodge", width = 0.75, color='black') +
  theme_classic() +
  labs(x='Tissue', y='Number of expressed HRs', fill='AGE') +
  scale_fill_brewer(palette = 'Blues') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave('output/HPA-single_cell/figure3/gtex-number-of-expressed-receptors-AGE.pdf', width = 12, height = 4)

#################################
# Number of expressed tissue for each HR（分局部和广泛表达类型）
df <- gtex_hr %>% group_by(Description, SMTS, Expression_universality) %>% summarise(TPM=mean(TPM)) # 计算HR在组织（SMTS）中的表达平均值
df <- df %>% group_by(Description, Expression_universality) %>% summarise(count = sum(log1p(TPM) > 0.1)) # 计算组织所表达的HR数目
mycoms <- list(c('selectively expressed HRs', 'ubiquitously expressed HRs'))
ggplot(df, aes(Expression_universality, count)) +
  geom_boxplot(aes(fill = Expression_universality), outlier.shape = NA) +
  #geom_jitter(width = 0.3) +
  theme_classic() +
  theme(legend.position = 'none') +
  stat_compare_means(method = 'wilcox.test', comparisons = mycoms) +
  labs(x=NULL, y='Number of tissues\nexpressing a certain HR')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
# ggsave('output/HPA-single_cell/figure3/gtex-number-of-expressed-tissues-universality.pdf', width = 1.5, height = 3)

#################################
# Number of expressed tissue for each HR（分核膜类型）
df <- gtex_hr %>% group_by(Description, SMTS, Receptor_subcellular) %>% summarise(TPM=mean(TPM)) # 计算HR在组织（SMTS）中的表达平均值
df <- df %>% group_by(Description, Receptor_subcellular) %>% summarise(count = sum(log1p(TPM) > 0.1)) # 计算组织所表达的HR数目
mycoms <- list(c('Membrane', 'Nucleus'))
ggplot(df, aes(Receptor_subcellular, count)) +
  geom_boxplot(aes(fill = Receptor_subcellular), outlier.shape = NA) +
  geom_jitter(width = 0.2, size=0.5) +
  theme_classic() +
  scale_fill_manual(values = c("#ED1B28","#237EB6"))+
  theme(legend.position = 'none') +
  stat_compare_means(method = 'wilcox.test', comparisons = mycoms) +
  labs(x=NULL, y='No. of expressed tissues')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure3/gtex-number-of-expressed-tissues-subcellular-1.pdf', width = 2, height = 4)

#################################
# 看具体是那些HR与Bladder和Adipose Tissue的衰老有关
# -----------------Bladder-----------------
age <- c('20-29','30-39','40-49','50-59','60-69','70-79')
df <- gtex_hr %>% group_by(Description, SMTS, AGE) %>% summarise(TPM=mean(TPM)) %>% ungroup() # 计算HR在组织（SMTS）中的表达平均值
df <- df %>% filter(SMTS == 'Bladder',log1p(TPM) > 0.1)
expressed_hr1 <- df %>% filter(AGE == '20-29') %>% pull(Description)
zz <- c()
for (i in 2:length(age)) {
  expressed_hr2 <- df %>% filter(AGE == age[i]) %>% pull(Description)
  print('-------------------------')
  print(age[i])
  diff1 <- setdiff(expressed_hr1, expressed_hr2) # 比上一个年龄段少了哪些HRs
  print(diff1)
  if(i != length(age)){
    zz <- c(zz, diff1)  
  }
  # diff2 <- setdiff(expressed_hr2, expressed_hr1) # 比上一个年龄段多了哪些HRs
  # print(diff2)
  expressed_hr1 <- expressed_hr2
}
# [1] "-------------------------"
# [1] "30-39"
# [1] "FSHR"   "HCRTR1" "LHCGR"  "MLNR"   "NPY4R" 
# [1] "-------------------------"
# [1] "40-49"
# [1] "AVPR1B" "GALR2"  "GLP1R"  "NR1I2"  "NR2E1"  "OPRD1"  "RXFP1" 
# [1] "-------------------------"
# [1] "50-59"
# [1] "FSHR"  "HNF4A" "RXFP4" "SCTR" 
# [1] "-------------------------"
# [1] "60-69"
# [1] "GALR2"  "GUCY2C" "HNF4G"  "MC4R"   "RXFP1"  "SSTR5" 

# -----------------Adipose Tissue-----------------
age <- c('20-29','30-39','40-49','50-59','60-69','70-79')
df <- gtex_hr %>% group_by(Description, SMTS, AGE) %>% summarise(TPM=mean(TPM)) %>% ungroup() # 计算HR在组织（SMTS）中的表达平均值
df <- df %>% filter(SMTS == 'Adipose Tissue',log1p(TPM) > 0.1)
expressed_hr1 <- df %>% filter(AGE == '20-29') %>% pull(Description)
for (i in 2:length(age)) {
  expressed_hr2 <- df %>% filter(AGE == age[i]) %>% pull(Description)
  print('-------------------------')
  print(age[i])
  # diff1 <- setdiff(expressed_hr1, expressed_hr2) # 比上一个年龄段少了哪些HRs
  # print(diff1)
  diff2 <- setdiff(expressed_hr2, expressed_hr1) # 比上一个年龄段多了哪些HRs
  print(diff2)
  expressed_hr1 <- expressed_hr2
}
# [1] "-------------------------"
# [1] "30-39"
# [1] "NR5A1" "SSTR5"
# [1] "-------------------------"
# [1] "40-49"
# [1] "CALCR" "ESRRG" "MC5R"  "NPY4R" "SSTR3"
# [1] "-------------------------"
# [1] "50-59"
# [1] "AGTR2" "CCKBR" "GLP1R" "NR5A1"
# [1] "-------------------------"
# [1] "60-69"
# [1] "GHRHR"
# [1] "-------------------------"
# [1] "70-79"
# [1] "MC5R"  "NR0B2" "NR5A1"

ggplot(df %>% filter(Description == 'GUCY2C'), aes(AGE, TPM)) +
  geom_bar(stat = 'identity') +
  theme_classic()


#################################
# HRs of differencet hormone classes distribution in tissues
gtex_hr <- gtex_hr %>% group_by(Description) %>% filter(TPM > (median(TPM) - 3*sd(TPM)), TPM < median(TPM) + 3*sd(TPM))
df <- gtex_hr %>% group_by(Description, SMTS, Hormone_chemical_classes) %>% summarise(TPM=mean(TPM)) # 计算HR在组织（SMTS）中的表达平均值
# df <- df %>% group_by(SMTS, Hormone_chemical_classes) %>% summarise(count = sum(log1p(TPM) > 0.1)) # 计算组织所表达的HR数目
# level <- unique(arrange(df %>% group_by(SMTS) %>% summarise(count = sum(count)), count)$SMTS)
# df$SMTS <- factor(df$SMTS, levels = level)

mycoms <- list(c('Amino acid-derived', 'Lipid-derived'), c('Lipid-derived', 'Peptide'))
ggplot(df, aes(Hormone_chemical_classes, log10(TPM+1), fill=Hormone_chemical_classes)) +
  facet_wrap(~SMTS, nrow = 3, strip.position = 'bottom') +
  geom_boxplot() +
  # theme_classic() +
  stat_summary(fun.y="mean", color='white', shape=3) +
  labs(x='Tissue', y='Log10 (TPM+1)', fill='Hormone classes') +
  ggpubr::stat_compare_means(method = 'wilcox.test', comparisons = mycoms, label = 'p.signif') +
  theme_classic() +
  theme(legend.position = 'top', axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank())
# ggsave('output/HPA-single_cell/figure3/gtex-hormone_class-2.pdf', width = 9.6, height = 10.6)


ggplot(df, aes(Hormone_chemical_classes, log10(TPM+1), fill=Hormone_chemical_classes)) +
  facet_wrap(~SMTS, nrow = 3, strip.position = 'bottom') +
  geom_violin() +
  # theme_classic() +
  stat_summary(fun.y="mean", color='white', shape=3) +
  stat_summary(fun.y="median", color='white', shape=20) +
  labs(x='Tissue', y='Log10 (TPM+1)', fill='Hormone classes') +
  ggpubr::stat_compare_means(method = 'wilcox.test', comparisons = mycoms, label = 'p.signif') +
  theme_classic() +
  theme(legend.position = 'top', axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank())
# ggsave('output/HPA-single_cell/figure3/gtex-hormone_class-violin.pdf', width = 10, height = 10.6)
