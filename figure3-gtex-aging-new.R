#################################
# Author: Xiangjie Zhao
# Date: 2024.07.09
# Title: Analysis of age-associated HRs
# using linear mixed model (LMM) or Pearson's cor test
#################################

rm(list = ls())
set.seed(1)
library(tidyverse)
source('./summarySE.R')

#################################
# showing the tissue expression of the age-related HRs.
gtex_hr <- read_tsv('../input/gtex/GTEx_HR_combined.tsv')
gtex_hr$AGE_raw <- gtex_hr$AGE
gtex_hr$AGE <- gtex_hr$AGE %>% recode(`20-29`=1, `30-39`=2, `40-49`=3, `50-59`=4, `60-69`=5, `70-79`=6)
# gtex_hr$TPM <- log2(gtex_hr$TPM)
tissue <- gtex_hr$SMTS%>% unique
gene <- gtex_hr$Description %>% unique

#################################
# # Pearson cor test
# df <- tibble(Gene = 'aaa', Tissue = 'bbb', `Pearson's Cor` = 0.3, `P value` = 0.05)
# for (i in tissue) {
#   for (j in gene) {
#     gtex_hr_ij <- gtex_hr %>% filter(SMTS == i, Description == j)
#     res <- cor.test(gtex_hr_ij$TPM, gtex_hr_ij$AGE, method = 'pearson')
#     df <- df %>% add_row(Gene=j, Tissue=i, `Pearson\'s Cor`=res$estimate, `P value`=res$p.value)
#   }
# }
# df <- df[-1,]
# df$P.adj <- p.adjust(df$`P value`, "BH")
# df <- df %>% arrange(Tissue, P.adj)
# write_tsv(df, '../output/HPA-single_cell/figure-aging-new/tissue-age-pearson-test-BH.tsv')
# df %>% filter(P.adj < 0.05) %>% write_tsv('../output/HPA-single_cell/figure-aging-new/tissue-age-pearson-test-BH-P.adj-0.05.tsv')
# df %>% filter(P.adj < 0.05) %>% filter(`Pearson's Cor` > 0.2 | `Pearson's Cor` < -0.2) %>% write_tsv('../output/HPA-single_cell/figure-aging-new/tissue-age-pearson-test-BH-P.adj-0.05-r0.2.tsv')


# # Pearson cor test (TPM cutoff-min)
# df <- tibble(Gene = 'aaa', Tissue = 'bbb', TPM_min = 0.7, `Pearson's Cor` = 0.3, `P value` = 0.05)
# for (i in tissue) {
#   for (j in gene) {
#     gtex_hr_ij <- gtex_hr %>% filter(SMTS == i, Description == j)
#     res <- cor.test(gtex_hr_ij$TPM, gtex_hr_ij$AGE, method = 'pearson')
#     tpm_min <- gtex_hr_ij %>% group_by(AGE_raw) %>% summarise(m = mean(TPM)) %>% pull(m) %>% min
#     df <- df %>% add_row(Gene=j, Tissue=i, TPM_min=tpm_min, `Pearson\'s Cor`=res$estimate, `P value`=res$p.value)
#   }
# }
# df <- df[-1,]
# df$P.adj <- p.adjust(df$`P value`, "BH")
# df <- df %>% arrange(Tissue, P.adj)
# write_tsv(df, '../output/HPA-single_cell/figure-aging-new/tissue-age-pearson-test-BH-withTPM.tsv')
# df %>% filter(P.adj < 0.05) %>% filter(`Pearson's Cor` > 0.2 | `Pearson's Cor` < -0.2) %>% filter(TPM_min > 1) %>% write_tsv('../output/HPA-single_cell/figure-aging-new/tissue-age-pearson-test-BH-P.adj-0.05-r0.2-TPMmin10.tsv')

############################
# # Pearson cor test (TPM cutoff-mean)
# df <- tibble(Gene = 'aaa', Tissue = 'bbb', tpm_mean = 0.7, `Pearson's Cor` = 0.3, `P value` = 0.05)
# for (i in tissue) {
#   for (j in gene) {
#     gtex_hr_ij <- gtex_hr %>% filter(SMTS == i, Description == j)
#     res <- cor.test(gtex_hr_ij$TPM, gtex_hr_ij$AGE, method = 'pearson')
#     tpm_mean <- gtex_hr_ij %>% group_by(AGE_raw) %>% summarise(m = mean(TPM)) %>% pull(m) %>% mean
#     df <- df %>% add_row(Gene=j, Tissue=i, tpm_mean=tpm_mean, `Pearson\'s Cor`=res$estimate, `P value`=res$p.value)
#   }
# }
# df <- df[-1,]
# df$P.adj <- p.adjust(df$`P value`, "BH")
# df <- df %>% arrange(Tissue, P.adj)
# # write_tsv(df, '../output/HPA-single_cell/figure-aging-new/tissue-age-pearson-test-BH-withTPM_mean.tsv')
# # df %>% filter(P.adj < 0.05) %>% filter(`Pearson's Cor` > 0.2 | `Pearson's Cor` < -0.2) %>% filter(tpm_mean > 1) %>% write_tsv('../output/HPA-single_cell/figure-aging-new/tissue-age-pearson-test-BH-P.adj-0.05-r0.2-TPMmean1.tsv')


#################################
# Rank tissue by number of significantly age-associated HRs
df <- readxl::read_excel('../output/HPA-single_cell/figure-aging-new/tissue-age-pearson-test-BH-P.adj-0.05-r0.2-TPMmean1.xlsx')
df_receptor <- readxl::read_excel('../input/hormone_receptor_list.xlsx')
df_receptor$Gene <- df_receptor$Gene_symbol
df <- left_join(df, df_receptor, by='Gene')
l <- df %>% group_by(Tissue) %>% count() %>% arrange(n) %>% pull(Tissue) # tissue ranks by # of age-related HRs

# df %>% write_tsv('../output/HPA-single_cell/figure-aging-new/tissue-age-pearson-test-BH-P.adj-0.05-r0.2-TPMmean1-new.tsv')

# group by HR types
df_count <- df %>% group_by(Tissue, Receptor_subcellular) %>% count() %>% arrange(n)
df_count$Tissue <- factor(df_count$Tissue, levels = l)
ggplot(df_count, aes(Tissue, n, fill=Receptor_subcellular)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  coord_flip() +
  labs(x='Tissue', y='No. of age-associated HRs', fill='HR types') +
  scale_fill_manual(values = c("#ED1B28","#237EB6")) +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('../output/HPA-single_cell/figure-aging-new/tissue_rank-hrType.pdf', width = 3.5, height = 4)

# group by hormone types
df_count <- df %>% group_by(Tissue, Hormone_chemical_classes) %>% count() %>% arrange(n)
df_count$Hormone_chemical_classes <- df_count$Hormone_chemical_classes %>% recode('NA'='Orphan')
df_count$Tissue <- factor(df_count$Tissue, levels = l)
ggplot(df_count, aes(Tissue, n, fill=Hormone_chemical_classes)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  coord_flip() +
  labs(x='Tissue', y='No. of age-associated HRs', fill='Hormone classes') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('../output/HPA-single_cell/figure-aging-new/tissue_rank-hormoneType.pdf', width = 4, height = 4)

################################
# show the top1 age-associated gene for each tissue using line plot with error bar
tissue <- df$Tissue %>% unique()
for (t in tissue) {
  g <- df %>% filter(Tissue == t) %>% head(1) %>% pull(Gene) # choose top 1 gene to show
  dd <- filter(gtex_hr, Description==g, SMTS==t)
  tgc <- summarySE(dd, measurevar='TPM', groupvars=c('AGE_raw'))
  ggplot(tgc, aes(x=AGE_raw, y=`TPM`)) +
    geom_errorbar(aes(ymin=`TPM`-se, ymax=`TPM`+se), color='black', width=0) +
    geom_line(group=1) +
    geom_point(size=1.5, shape=21, fill='white', stroke=0.7) + # 21 is filled circle
    labs(x='Age', y=paste0(g,' (TPM)'), title = t) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.text = element_text(color="black"),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5))
  ggsave(paste0('../output/HPA-single_cell/figure-aging-new/new-TPM-',t,'-',g,'.pdf'), width = 1.5, height = 2)
}


# 
# t <- 'Fallopian Tube'
# g <- 'HCRTR1'
# dd <- filter(gtex_hr, Description==g, SMTS==t)
# tgc <- summarySE(dd, measurevar='TPM', groupvars=c('AGE_raw'))
# ggplot(tgc, aes(x=AGE_raw, y=`TPM`)) +
#   geom_errorbar(aes(ymin=`TPM`-se, ymax=`TPM`+se), color='black', width=0) +
#   geom_line(group=1) +
#   geom_point(size=1.5, shape=21, fill='white', stroke=0.7) + # 21 is filled circle
#   labs(x='Age', y=paste0(g,' (TPM)'), title = t) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         axis.text = element_text(color="black"),
#         axis.ticks = element_line(color = "black"),
#         plot.title = element_text(hjust = 0.5))
# ggsave(paste0('../output/HPA-single_cell/figure-aging-new/TPM-',t,'-',g,'.pdf'), width = 1.5, height = 2)


#################################
# 分析年龄相关HR的普适性：HR按照在组织中显著出现的次数排序（min 0, max 30）
df <- readxl::read_excel('../output/HPA-single_cell/figure-aging-new/tissue-age-pearson-test-BH-P.adj-0.05-r0.2-TPMmean1.xlsx')
df_gene_count <- df %>% group_by(Gene) %>% count() %>% arrange(n) %>% tail(10) # show 10 HRs
ggplot(df_gene_count, aes(factor(Gene, levels = rev(df_gene_count$Gene)), n)) +
  geom_bar(stat = 'identity', fill='#8AA0C9', width=0.5) +
  theme_classic() +
  # coord_flip() +
  labs(x='Gene', y='No. of tissues') +
  scale_y_continuous(breaks = c(0, 3, 6, 9)) +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave('../output/HPA-single_cell/figure-aging-new/HR-rank-by-No.tissue.pdf', width = 5, height = 2)


#################################
# 组织间的相关性。分析组织间share的显著age相关HR。用matrix热图展示。（或者上三角）。X, Y分别都是组织类型，中间对应的数值是overlap的HR个数
df <- readxl::read_excel('../output/HPA-single_cell/figure-aging-new/tissue-age-pearson-test-BH-P.adj-0.05-r0.2-TPMmean1.xlsx')
tissue <- df$Tissue %>% unique
df_long <- tibble(Tissue1 = 'aaa', Tissue2 = 'bbb', `No. of shared HRs` = 7)
for (t1 in tissue) {
  for (t2 in tissue) {
    hr_t1 <- df %>% filter(Tissue == t1) %>% pull(Gene)
    hr_t2 <- df %>% filter(Tissue == t2) %>% pull(Gene)
    n <- length(intersect(hr_t1, hr_t2))
    df_long <- df_long %>% add_row(Tissue1=t1, Tissue2=t2, `No. of shared HRs`=n)
  }
}
df_long <- df_long[-1,]

# keep bottom triangle of the matrix
df_long <- df_long %>% rowwise() %>%
  mutate(pair = sort(c(Tissue1, Tissue2)) %>% paste(collapse = ",")) %>%
  group_by(pair) %>%
  distinct(pair, .keep_all = T) %>%
  filter(!pair %in% str_c(tissue,',',tissue)) # drop pairs like "tissue A-tissue A"
# head(df_long)

# plotting
ggplot(df_long, aes(x = Tissue1, y = Tissue2, fill = `No. of shared HRs`)) +
  geom_raster() +
  geom_text(aes(label = `No. of shared HRs`), size = 2.5) +
  scale_fill_distiller(palette = "Spectral") +
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust=0),
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"))
ggsave('../output/HPA-single_cell/figure-aging-new/triangle-heatmap-No. of shared HRs.pdf', width = 7.5, height = 5.4)


#################################
# 组织间的相关性。分析组织间share的显著age相关HR。用matrix热图展示。（或者上三角）。X, Y分别都是组织类型，中间对应的数值是overlap的HR个数
# 中间数值改为overlap的HR个数处以normalized factor。normalized factor是对应两个组织age-associated HR数目的均值。
# 观察上面的图，发现热点区域都是本身age-associated HR数目多的组织，这样标准化可以减少组织本身拥有age-associated HR数目的影响。
df <- readxl::read_excel('../output/HPA-single_cell/figure-aging-new/tissue-age-pearson-test-BH-P.adj-0.05-r0.2.xlsx')
tissue <- df$Tissue %>% unique
df_long <- tibble(Tissue1 = 'aaa', Tissue2 = 'bbb', `No. of shared HRs` = 7)
for (t1 in tissue) {
  for (t2 in tissue) {
    hr_t1 <- df %>% filter(Tissue == t1) %>% pull(Gene)
    hr_t2 <- df %>% filter(Tissue == t2) %>% pull(Gene)
    n <- length(intersect(hr_t1, hr_t2))
    df_long <- df_long %>% add_row(Tissue1=t1, Tissue2=t2, `No. of shared HRs`=n)
  }
}
df_long <- df_long[-1,]

v1 <- df %>% group_by(Tissue) %>% count() %>% pull(n) %>% rep(30)
v2 <- df %>% group_by(Tissue) %>% count() %>% pull(n) %>% rep(each=30)
v <- (v1 + v2) / 2

df_long$`No. of shared HRs` <- (df_long$`No. of shared HRs`  / v) %>% round(1)

# keep bottom triangle of the matrix
df_long <- df_long %>% rowwise() %>%
  mutate(pair = sort(c(Tissue1, Tissue2)) %>% paste(collapse = ",")) %>%
  group_by(pair) %>%
  distinct(pair, .keep_all = T) %>%
  filter(!pair %in% str_c(tissue,',',tissue)) # drop pairs like "tissue A-tissue A"
# head(df_long)

# plotting
ggplot(df_long, aes(x = Tissue1, y = Tissue2, fill = `No. of shared HRs`)) +
  geom_raster() +
  geom_text(aes(label = `No. of shared HRs`), size = 2.5) +
  scale_fill_distiller(palette = "Spectral") +
  theme_minimal() +
  scale_x_discrete(position = "top") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust=0),
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"))
# ggsave('../output/HPA-single_cell/figure-aging-new/triangle-heatmap-No. of shared HRs-normalized.pdf', width = 9, height = 7)


