#################################
# Author: Xiangjie Zhao
# Date: 2023.08.07
# Title: Analysis of hormone receptor gene list (SEX)
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
sex_tissue <- c('Prostate','Testis','Cervix Uteri','Fallopian Tube','Ovary','Uterus','Vagina')

#################################
gtex_hr <- read_tsv('input/gtex/GTEx_HR_combined.tsv')
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
df_receptor$Hormone_chemical_classes <- recode(df_receptor$Hormone_chemical_classes, 'NA'='Orphan')
df_receptor$Description <- df_receptor$Gene_symbol
gtex_hr <- left_join(gtex_hr, df_receptor, by='Description')

# ################################
# # 设置阈值，基因在组织中的平均TPM低于此值视为不表达
# df <- gtex_hr %>% group_by(Description, SMTS) %>% summarise(TPM=mean(TPM))
# # 这里即log1p(TPM) < 0.1视为不表达
# ggplot(df, aes(log1p(TPM)))+
#   geom_density()+
#   geom_vline(xintercept=0.1)+
#   labs(y='Distribution density')+
#   theme_classic()
# # ggsave('output/HPA-single_cell/figure3/gtex-cutoff.pdf', width = 2.5, height = 2)
# 
# #################################
# # 30种组织表达HR数目的排名（分HR核膜类型）
# df <- gtex_hr %>% group_by(Description, SMTS, Receptor_subcellular) %>% summarise(TPM=mean(TPM)) # 计算HR在组织（SMTS）中的表达平均值
# df <- df %>% group_by(SMTS, Receptor_subcellular) %>% summarise(count = sum(log1p(TPM) > 0.1)) # 计算组织所表达的HR数目
# level <- unique(arrange(df %>% group_by(SMTS) %>% summarise(count = sum(count)), count)$SMTS)
# df$SMTS <- factor(df$SMTS, levels = level)
# ggplot(df, aes(SMTS, count, fill=Receptor_subcellular)) +
#   geom_bar(stat = 'identity') +
#   theme_classic() +
#   # coord_flip() +
#   labs(x='Tissue', y='No. of HRs', fill='HR types') +
#   scale_y_continuous(breaks = seq(0,140,10)) +
#   scale_fill_manual(values = c("#ED1B28","#237EB6")) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# # ggsave('output/HPA-single_cell/figure3/gtex-number-of-expressed-receptors-HRtype-1.pdf', width = 7, height = 5)
# 
# #################################
# # 30种组织表达HR数目的排名（分性别）
# df <- gtex_hr %>% group_by(Description, SMTS, SEX) %>% summarise(TPM=mean(TPM)) # 计算HR在组织（SMTS）中的表达平均值
# df <- df %>% group_by(SMTS, SEX) %>% summarise(count = sum(log1p(TPM) > 0.1)) # 计算组织所表达的HR数目
# level <- unique(arrange(df %>% group_by(SMTS) %>% summarise(count = sum(count)), count)$SMTS)
# df$SMTS <- factor(df$SMTS, levels = level)
# ggplot(df, aes(SMTS, count, fill=SEX)) +
#   geom_bar(stat = 'identity', position = "dodge", width = 0.7) +
#   theme_classic() +
#   labs(x='Tissue', y='Number of expressed HRs', fill='SEX') +
#   scale_fill_manual(values = c("#3BAF51","#9C4FA1")) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# # ggsave('output/HPA-single_cell/figure-sex/gtex-number-of-expressed-receptors-SEX.pdf', width = 8, height = 4)
# 
# # 男性相比女性各种组织所表达HR数目变化的百分比（负值代表减少）
# df_female <- df %>% filter(SEX == 'Female', SMTS %in% common_tissue)
# df_male <- df %>% filter(SEX == 'Male', SMTS %in% common_tissue)
# df_female$change_proportion <- (df_male$count - df_female$count) / df_female$count * 100
# df_female$change_direction <- to_vec(for(i in df_female$change_proportion) if(i<0) 'Down' else 'Up')
# ggplot(df_female, aes(SMTS, change_proportion, fill=change_direction)) +
#   geom_bar(stat = 'identity', color='black', linewidth=0.3) +
#   coord_flip() +
#   theme_classic() +
#   labs(x='Tissue', y='Changed percentage', fill='Change direction') +
#   scale_fill_manual(values = c('#237EB5','#EE1B27'))
# # ggsave('output/HPA-single_cell/figure-sex/gtex-number-of-expressed-receptors-SEX-changePercentageMale2Female.pdf', width = 4, height = 4)
# 
# # 具体查看Adipose Tissue中男性比女性少表达的HR是哪几个。
# hr_expressed <- gtex_hr %>% group_by(Description, SMTS, SEX) %>% summarise(TPM=mean(TPM)) %>% filter(log1p(TPM) > 0.1) %>% filter(SMTS == 'Thyroid')
# hr_expressed_male <- hr_expressed %>% filter(SEX == 'Male') %>% pull(Description)
# hr_expressed_female <- hr_expressed %>% filter(SEX == 'Female') %>% pull(Description)
# setdiff(hr_expressed_male, hr_expressed_female)
# setdiff(hr_expressed_female, hr_expressed_male)

#################################
# # Wilcox test for TPM (male v.s. female, non-sex tissue)
# gtex_hr_common <- gtex_hr %>% filter(SMTS %in% common_tissue) # only test non-sex tissue
# tissue <- gtex_hr_common$SMTS%>% unique
# gene <- gtex_hr_common$Description %>% unique
# 
# df <- tibble(Gene = 'aaa', Tissue = 'bbb', `Log2 foldchange` = 0.3, `P value` = 0.05)
# for (i in tissue) {
#   for (j in gene) {
#     gtex_hr_ij <- gtex_hr_common %>% filter(SMTS == i, Description == j)
#     gtex_hr_ij_male <- filter(gtex_hr_ij, SEX == 'Male')
#     gtex_hr_ij_female <- filter(gtex_hr_ij, SEX == 'Female')
#     p_value <- wilcox.test(gtex_hr_ij_male$TPM, gtex_hr_ij_female$TPM)$p.value
#     lgfc <- log2(mean(gtex_hr_ij_male$TPM) / mean(gtex_hr_ij_female$TPM))
#     df <- df %>% add_row(Gene=j, Tissue=i, `Log2 foldchange`=lgfc, `P value`=p_value)
#   }
# }
# df <- df[-1,]
# df$P.adj <- p.adjust(df$`P value`, "BH")
# df <- df %>% arrange(Tissue, P.adj)
# write_tsv(df, 'output/HPA-single_cell/figure-sex/tissue-sex-wilcox-test-BH.tsv')
# df %>% filter(P.adj < 0.05) %>% filter((`Log2 foldchange` > log2(1.5)) | (`Log2 foldchange` < log2(1/1.5))) %>% write_tsv('output/HPA-single_cell/figure-sex/tissue-sex-wilcox-test-BH-P.adj-0.05-fc-1.5.tsv')

# 
# # Wilcox test for TPM (male v.s. female, non-sex tissue)-平均TPM
# gtex_hr_common <- gtex_hr %>% filter(SMTS %in% common_tissue) # only test non-sex tissue
# tissue <- gtex_hr_common$SMTS%>% unique
# gene <- gtex_hr_common$Description %>% unique
# 
# df <- tibble(Gene = 'aaa', Tissue = 'bbb', TPM_mean = 0.7,`Log2 foldchange` = 0.3, `P value` = 0.05)
# for (i in tissue) {
#   for (j in gene) {
#     gtex_hr_ij <- gtex_hr_common %>% filter(SMTS == i, Description == j)
#     gtex_hr_ij_male <- filter(gtex_hr_ij, SEX == 'Male')
#     gtex_hr_ij_female <- filter(gtex_hr_ij, SEX == 'Female')
#     p_value <- wilcox.test(gtex_hr_ij_male$TPM, gtex_hr_ij_female$TPM)$p.value
#     lgfc <- log2(mean(gtex_hr_ij_male$TPM) / mean(gtex_hr_ij_female$TPM))
#     df <- df %>% add_row(Gene=j, Tissue=i, TPM_mean=(mean(gtex_hr_ij_male$TPM)+mean(gtex_hr_ij_female$TPM))/2, `Log2 foldchange`=lgfc, `P value`=p_value)
#   }
# }
# df <- df[-1,]
# df$P.adj <- p.adjust(df$`P value`, "BH")
# df <- df %>% arrange(Tissue, P.adj)
# write_tsv(df, 'output/HPA-single_cell/figure-sex/tissue-sex-wilcox-test-BH-TPMmean.tsv')
# # df %>% filter(P.adj < 0.05) %>% filter((`Log2 foldchange` > log2(1.5)) | (`Log2 foldchange` < log2(1/1.5))) %>% filter(TPM_mean > 1) %>% write_tsv('output/HPA-single_cell/figure-sex/tissue-sex-wilcox-test-BH-P.adj-0.05-fc-1.5-TPMmean1.tsv')


#################################
# 看下组织的性别差异HR的总体情况
df <- readxl::read_excel('output/HPA-single_cell/figure-sex/tissue-sex-wilcox-test-BH-P.adj-0.05-fc-1.5-TPMmean1.xlsx')
library(ggrepel)
library(randomcoloR)
set.seed(3)
c <- distinctColorPalette(14)
ggplot(df, aes(`Log2 foldchange`, -log10(`P value`), fill=Tissue, color=Tissue, label=Gene)) +
  geom_text_repel(max.overlaps = 100, size=2.2) +
  geom_point(shape=21, size=1.5, stroke=0.1, color='black') +
  scale_fill_manual(values = c) +
  scale_color_manual(values = c) +
  theme_bw() +
  labs(x='Log2 FoldChange', y='-log10 (P value)')+
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure-sex/volcano-sexDiffHR-new.pdf', width = 5, height = 8)

#################################
# Rank tissue by number of significantly sex-associated HRs
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
df_receptor$Gene <- df_receptor$Gene_symbol
df <- left_join(df, df_receptor, by='Gene')
l <- df %>% group_by(Tissue) %>% count() %>% arrange(n) %>% pull(Tissue) # tissue ranks by # of sex-related HRs

# df %>% write_tsv('output/HPA-single_cell/figure-sex/tissue-sex-wilcox-test-BH-P.adj-0.05-fc-1.5-TPMmean1-new.tsv')

# no grouping
df_count <- df %>% group_by(Tissue) %>% count() %>% arrange(n)
df_count$Tissue <- factor(df_count$Tissue, levels = l)
ggplot(df_count, aes(Tissue, n)) +
  geom_bar(stat = 'identity', color='black', fill='white') +
  theme_classic() +
  coord_flip() +
  labs(x='Tissue', y='No. of sex-associated HRs') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure-sex/tissue_rank.pdf', width = 3, height = 2.5)

# group by HR types
df_count <- df %>% group_by(Tissue, Receptor_subcellular) %>% count() %>% arrange(n)
df_count$Tissue <- factor(df_count$Tissue, levels = l)
ggplot(df_count, aes(Tissue, n, fill=Receptor_subcellular)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  coord_flip() +
  labs(x='Tissue', y='No. of sex-associated HRs', fill='HR types') +
  scale_fill_manual(values = c("#ED1B28","#237EB6")) +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
ggsave('output/HPA-single_cell/figure-sex/tissue_rank-hrType.pdf', width = 4, height = 2.5)

# group by hormone types
df_count <- df %>% group_by(Tissue, Hormone_chemical_classes) %>% count() %>% arrange(n)
df_count$Hormone_chemical_classes <- df_count$Hormone_chemical_classes %>% recode('NA'='Orphan')
df_count$Tissue <- factor(df_count$Tissue, levels = l)
ggplot(df_count, aes(Tissue, n, fill=Hormone_chemical_classes)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  coord_flip() +
  labs(x='Tissue', y='No. of sex-associated HRs', fill='Hormone classes') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
ggsave('output/HPA-single_cell/figure-sex/tissue_rank-hormoneType.pdf', width = 4.5, height = 2.5)

# group by sex
df_recode <- df
df_recode$sex <- if_else(df$`Log2 foldchange` > 0, 'Male-biased', 'Female-biased')
df_count <- df_recode %>% group_by(Tissue, sex) %>% count() %>% arrange(n)
df_count$sex <- df_count$sex %>% recode('NA'='Orphan')
df_count$Tissue <- factor(df_count$Tissue, levels = l)
ggplot(df_count, aes(Tissue, n, fill=sex)) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  coord_flip() +
  labs(x='Tissue', y='No. of sex-associated HRs', fill='Expression preference') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure-sex/tissue_rank-sexBiased.pdf', width = 4.5, height = 2.5)

# group by sex-1
df_recode <- df
df_recode$sex <- if_else(df$`Log2 foldchange` > 0, 'Male-biased', 'Female-biased')
df_count <- df_recode %>% group_by(Tissue, sex) %>% count() %>% arrange(n)
df_count$sex <- df_count$sex %>% recode('NA'='Orphan')
df_count$Tissue <- factor(df_count$Tissue, levels = l)
ggplot(df_count, aes(Tissue, n, fill=sex)) +
  geom_bar(stat = 'identity', position = 'stack') +
  theme_classic() +
  coord_flip() +
  labs(x='Tissue', y='No. of sex-associated HRs', fill='Expression preference') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
ggsave('output/HPA-single_cell/figure-sex/tissue_rank-sexBiased-1.pdf', width = 4.6, height = 2.5)


#################################
# 分析性别表达有差异HR的普适性：HR按照在组织中显著出现的次数排序
df_gene_count <- df %>% group_by(Gene) %>% count() %>% arrange(n) %>% filter(n > 1) # show n >= 2 HRs
ggplot(df_gene_count, aes(factor(Gene, levels = rev(df_gene_count$Gene)), n)) +
  geom_bar(stat = 'identity', fill='#8AA0C9') +
  theme_classic() +
  # coord_flip() +
  labs(x='Gene', y='No. of tissues') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('output/HPA-single_cell/figure-sex/HR-rank-by-No.tissue-new.pdf', width = 3, height = 1.8)


###########################################################
# 画一下某种组织中这些基因的表达量图，看下男女表达差异。
# show the top1 (arrange by foldchange) sex-associated HR for each tissue using boxplot
df <- readxl::read_excel('output/HPA-single_cell/figure-sex/tissue-sex-wilcox-test-BH-P.adj-0.05-fc-1.5-TPMmean1.xlsx')
tissue <- df$Tissue %>% unique()
for (t in tissue) {
  hr_stats <- df %>% filter(Tissue == t) %>% arrange(abs(`Log2 foldchange`)) %>% tail(1) # choose top 1 gene to show
  hr_tpm <- filter(gtex_hr, Description==hr_stats$Gene, SMTS==t)
  hr_tpm <- hr_tpm %>% group_by(SEX) %>% filter(TPM < (mean(TPM) + 2*sd(TPM)), TPM > (mean(TPM) - 2*sd(TPM)))
  ggplot(hr_tpm, aes(SEX, TPM, fill=SEX)) +
    geom_boxplot() +
    scale_fill_manual(values = c("#3BAF51","#9C4FA1")) +
    theme_classic() +
    labs(x='', y=paste0(hr_stats$Gene,' (TPM)'), title = t,
         subtitle = paste0('p = ', signif(hr_stats$`P.adj`, 2))) +
    theme(legend.position = 'none',
          axis.text = element_text(color="black"),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5))
  # ggsave(paste0('output/HPA-single_cell/figure-sex/new-TPM-',t,'-',hr_stats$Gene,'.pdf'), width = 1.5, height = 3.5)
}

# CRHR1
t <- 'Pituitary'
hr_stats <- df %>% filter(Tissue == t) %>% filter(Gene == 'CRHR1')
hr_tpm <- filter(gtex_hr, Description==hr_stats$Gene, SMTS==t)
hr_tpm <- hr_tpm %>% group_by(SEX) %>% filter(TPM < (mean(TPM) + 2*sd(TPM)), TPM > (mean(TPM) - 2*sd(TPM)))
ggplot(hr_tpm, aes(SEX, TPM, fill=SEX)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#3BAF51","#9C4FA1")) +
  theme_classic() +
  labs(x='', y=paste0(hr_stats$Gene,' (TPM)'), title = t,
       subtitle = paste0('p = ', signif(hr_stats$`P.adj`, 2))) +
  theme(legend.position = 'none',
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5))
ggsave(paste0('output/HPA-single_cell/figure-sex/new-TPM-',t,'-',hr_stats$Gene,'.pdf'), width = 1.5, height = 3.5)




###########################################################
# 性激素受体的表达差异
g <- 'NR1I3'
hr_tpm <- filter(gtex_hr, Description==g)
hr_tpm <- hr_tpm %>% group_by(SEX, SMTS) %>% filter(TPM < (mean(TPM) + 2*sd(TPM)), TPM > (mean(TPM) - 2*sd(TPM)))
ggplot(hr_tpm, aes(SEX, TPM, fill=SEX)) +
  geom_boxplot(width=0.6) +
  facet_wrap(~SMTS) +
  scale_fill_manual(values = c("#3BAF51","#9C4FA1")) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure-sex/sex_HR-TPM-NR1I3.pdf', width = 6.5, height = 5)

#------------------------------
# example
g <- 'GPER1'
hr_tpm <- filter(gtex_hr, Description==g, SMTS=='Breast')
hr_tpm <- hr_tpm %>% group_by(SEX, SMTS) %>% filter(TPM < (mean(TPM) + 2*sd(TPM)), TPM > (mean(TPM) - 2*sd(TPM)))
ggplot(hr_tpm, aes(SEX, TPM, fill=SEX)) +
  geom_boxplot(width=0.6) +
  scale_fill_manual(values = c("#3BAF51","#9C4FA1")) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"))

#------------------------------
# 在性别器官中偏好表达的HR的分析

# 检验性器官VS非性器官表达差异高的基因
# 男性器官：male_tissue <- c('Prostate', 'Testis')
# 女性器官：female_tissue <- c('Cervix Uteri', 'Fallopian Tube', 'Ovary', 'Uterus', 'Vagina')
# 
# tissue <- c('Prostate', 'Testis', 'Cervix Uteri', 'Fallopian Tube', 'Ovary', 'Uterus', 'Vagina')
# gene <- gtex_hr$Description %>% unique
# df <- tibble(Gene = 'aaa', Tissue = 'bbb', TPM_sexTissue = 0.7,`Log2 foldchange` = 0.3, `P value` = 0.05)
# for (i in tissue) {
#   for (j in gene) {
#     gtex_hr_ij <- gtex_hr %>% filter(Description == j)
#     gtex_hr_ij_sex_tissue <- filter(gtex_hr_ij, SMTS == i)
#     gtex_hr_ij_nonsex_tissue <- filter(gtex_hr_ij, SEX != i)
#     p_value <- wilcox.test(gtex_hr_ij_sex_tissue$TPM, gtex_hr_ij_nonsex_tissue$TPM, alternative='greater')$p.value
#     lgfc <- log2(mean(gtex_hr_ij_sex_tissue$TPM) / mean(gtex_hr_ij_nonsex_tissue$TPM))
#     df <- df %>% add_row(Gene=j, Tissue=i, TPM_sexTissue=mean(gtex_hr_ij_sex_tissue$TPM), `Log2 foldchange`=lgfc, `P value`=p_value)
#   }
# }
# df <- df[-1,]
# df$P.adj <- p.adjust(df$`P value`, "BH")
# df <- df %>% arrange(Tissue, P.adj)
# # write_tsv(df, 'output/HPA-single_cell/figure-sex/sex2nonsex-wilcox-test-BH.tsv')
# # 选择其中高显著、高foldchange和高表达的HR
# # df %>% filter(P.adj < 0.05) %>% filter(`Log2 foldchange` > log2(2)) %>% write_tsv('output/HPA-single_cell/figure-sex/sex2nonsex-wilcox-test-BH-P.adj-0.05-fc-2.tsv')
# # df %>% filter(P.adj < 0.05) %>% filter(`Log2 foldchange` > log2(2)) %>% filter(log10(TPM_sexTissue) > 1) %>% write_tsv('output/HPA-single_cell/figure-sex/sex2nonsex-wilcox-test-BH-P.adj-0.05-fc-2-TPM1.tsv')


#-------
# tissue <- c('Prostate', 'Testis', 'Cervix Uteri', 'Fallopian Tube', 'Ovary', 'Uterus', 'Vagina')
# gene <- gtex_hr$Description %>% unique
# df <- tibble(Gene = 'aaa', Tissue = 'bbb', TPM_mean = 0.7,`Log2 foldchange` = 0.3, `P value` = 0.05)
# for (i in tissue) {
#   for (j in gene) {
#     gtex_hr_ij <- gtex_hr %>% filter(Description == j)
#     gtex_hr_ij_sex_tissue <- filter(gtex_hr_ij, SMTS == i)
#     gtex_hr_ij_nonsex_tissue <- filter(gtex_hr_ij, SEX != i)
#     p_value <- wilcox.test(gtex_hr_ij_sex_tissue$TPM, gtex_hr_ij_nonsex_tissue$TPM, alternative='greater')$p.value
#     lgfc <- log2(mean(gtex_hr_ij_sex_tissue$TPM) / mean(gtex_hr_ij_nonsex_tissue$TPM))
#     df <- df %>% add_row(Gene=j, Tissue=i, TPM_mean=(mean(gtex_hr_ij_nonsex_tissue$TPM)+mean(gtex_hr_ij_sex_tissue$TPM)) / 2, `Log2 foldchange`=lgfc, `P value`=p_value)
#   }
# }
# df <- df[-1,]
# df$P.adj <- p.adjust(df$`P value`, "BH")
# df <- df %>% arrange(Tissue, P.adj)
# write_tsv(df, 'output/HPA-single_cell/figure-sex/sex2nonsex-wilcox-test-BH-TPM_mean.tsv')
# # 选择其中高显著、高foldchange和高表达的HR
# df %>% filter(P.adj < 0.05) %>% filter(`Log2 foldchange` > log2(2)) %>% filter(TPM_mean > 3) %>% write_tsv('output/HPA-single_cell/figure-sex/sex2nonsex-wilcox-test-BH-P.adj-0.05-fc-2-TPM_mean3.tsv')


#-------

# # 设置cut-off大于多少就算性别相关HR。放在图1前面
# # 设置阈值，基因在组织中的平均TPM高于此值视为高表达
# zz <- gtex_hr %>% group_by(Description, SMTS) %>% summarise(TPM=mean(TPM))
# # 这里即log10(TPM) > 1视高表达
# ggplot(zz, aes(log10(TPM)))+
#   geom_density()+
#   geom_vline(xintercept=1)+
#   labs(y='Distribution density')+
#   theme_classic()
# ggsave('output/HPA-single_cell/figure-sex/gtex-cutoff-sex.pdf', width = 2.5, height = 2)

#------------------
df <- readxl::read_excel('output/HPA-single_cell/figure-sex/sex2nonsex-wilcox-test-BH-P.adj-0.05-fc-2-TPM_mean3.xlsx')
df <- df %>% arrange(Tissue, desc(`Log2 foldchange`))
df %>% group_by(Tissue) %>% top_n(2, `Log2 foldchange`)
df %>% group_by(Tissue) %>% top_n(1, `Log2 foldchange`)

df_plot <- df %>% group_by(Gene) %>% count() %>% arrange(desc(n))
# df_plot <- df %>% group_by(Gene) %>% count() %>% arrange(n) %>% filter(n > 1) # show n >= 2 HRs
ggplot(df_plot, aes(factor(Gene, levels = rev(df_plot$Gene)), n)) +
  geom_bar(stat = 'identity', fill='#8AA0C9') +
  theme_classic() +
  coord_flip() +
  labs(x='Gene', y='No. of tissues') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure-sex/sex organ-specific-count.pdf', width = 3.5, height = 7)


#------------------
# example
g <- 'MC1R'
hr_tpm <- filter(gtex_hr, Description==g)
hr_tpm <- hr_tpm %>% group_by(SEX, SMTS) %>% filter(TPM < (mean(TPM) + 2*sd(TPM)), TPM > (mean(TPM) - 2*sd(TPM)))
ggplot(hr_tpm, aes(SEX, TPM, fill=SEX)) +
  geom_boxplot(width=0.6) +
  facet_wrap(~SMTS) +
  scale_fill_manual(values = c("#3BAF51","#9C4FA1")) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"))
ggsave('output/HPA-single_cell/figure-sex/sex_HR-TPM-MC1R.pdf', width = 6.5, height = 5)
