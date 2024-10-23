#################################
# Author: Xiangjie Zhao
# Date: 2022.12.24
# Title: Analysis of human single-cell data (HPA single-cell data)
#################################
rm(list = ls())
set.seed(1)

library(tidyverse)
library(ComplexHeatmap)
library(tidyHeatmap)
library(ggpubr)
library(comprehenr)
library(RColorBrewer)
library(ggwordcloud)
setwd('~/project/hormone_receptor_scRNA-2022-12-11/')

# 74 colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#################################
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
vec_receptor <- df_receptor$Gene_symbol

df_celltype <- read.table('input/HPA-single_cell/single-cell-type/rna_single_cell_type.tsv',
                          sep = '\t', header = T)
df_celltype <- filter(df_celltype, Gene.name %in% vec_receptor)
df_celltype$Cell.type <- df_celltype$Cell.type %>% str_to_title


#################################
# mean TPM of membrane HR and nuclear HR
df_receptor <- df_receptor[str_order(df_receptor$Gene_symbol),]
df_mean_TPM <- df_celltype
df_mean_TPM <- df_mean_TPM[str_order(df_mean_TPM$Gene.name),]
df_mean_TPM$Receptor_subcellular <- rep(df_receptor$Receptor_subcellular, each=length(unique(df_mean_TPM$Cell.type)))
df_mean_TPM <- df_mean_TPM %>% group_by(Gene.name, Receptor_subcellular) %>% summarise(TPM=mean(nTPM))

mycoms <- list(c('Membrane', 'Nucleus'))
ggplot(df_mean_TPM, aes(Receptor_subcellular, TPM)) +
  geom_boxplot(aes(fill = Receptor_subcellular), outlier.shape = 1) +
  #geom_jitter(width = 0.3) +
  theme_classic() +
  scale_fill_manual(values = c("#ED1B28","#237EB6"))+
  theme(legend.position = 'none') +
  stat_compare_means(method = 'wilcox.test', comparisons = mycoms) +
  labs(x=NULL, y='TPM')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure2/mean-TPM-of-membrane-HR-and-nuclear-HR.pdf', width = 1.5, height = 2.5)


#################################
# Number of expressed hormone receptors in each cell types
df_count_1 <- df_celltype %>% group_by(Cell.type) %>% summarize(count = sum(nTPM > 0)) %>% arrange(count)
df_count_1$Cell.type <- factor(df_count_1$Cell.type, levels = df_count_1$Cell.type)

df_receptor <- df_receptor[str_order(df_receptor$Gene_symbol),]
df_count <- df_celltype
df_count <- df_count[str_order(df_count$Gene.name),]
df_count <- df_count[str_order(df_count$Cell.type),]
df_count$Receptor_subcellular <- rep(df_receptor$Receptor_subcellular, length(unique(df_celltype$Cell.type)))
df_count <- df_count %>% group_by(Cell.type, Receptor_subcellular) %>% summarize(count = sum(nTPM > 0)) %>% arrange(Cell.type, Receptor_subcellular)
df_count$Cell.type <- factor(df_count$Cell.type, levels = levels(df_count_1$Cell.type))

ggplot(df_count, aes(Cell.type, count, fill=Receptor_subcellular)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  coord_flip() +
  labs(x='Cell type', y='Number of expressed HRs', fill='HR types') +
  scale_fill_manual(values = c("#ED1B28","#237EB6"))

# ggsave('output/HPA-single_cell/figure2/number-of-expressed-receptors-cellType.pdf', width = 6, height = 10.5)


#################################
# Number of expressed hormone receptors in each cell type groups
df_cell_cluster_description <- read.table('input/HPA-single_cell/single-cell-type-tissue-cluster/rna_single_cell_cluster_description.tsv', sep = '\t', header = T)
df_cell_cluster_description$Cell.type <- df_cell_cluster_description$Cell.type %>% str_to_title
df_cell_cluster_description_group <- df_cell_cluster_description %>% select(Cell.type, Cell.type.group) %>% unique() %>% as.data.frame
df_celltype_group <- left_join(df_celltype, df_cell_cluster_description_group, by='Cell.type')

df_receptor$Gene.name <- df_receptor$Gene_symbol
df_celltype_group <- left_join(df_celltype_group, df_receptor, by='Gene.name')
df_count <- df_celltype_group %>% group_by(Cell.type.group, Gene.name, Receptor_subcellular) %>% summarize(nTPM = mean(nTPM)) %>% group_by(Cell.type.group, Receptor_subcellular) %>% summarize(count = sum(nTPM > 0)) %>% arrange(Cell.type.group, Receptor_subcellular)
level <- unique(arrange(df_count %>% group_by(Cell.type.group) %>% summarise(count = sum(count)), count)$Cell.type.group)
df_count$Cell.type.group <- factor(df_count$Cell.type.group, levels = level)
ggplot(df_count, aes(Cell.type.group, count, fill=Receptor_subcellular)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  coord_flip() +
  labs(x='Cell type group', y='Number of expressed HRs', fill='HR types') +
  scale_fill_manual(values = c("#ED1B28","#237EB6"))
# ggsave('output/HPA-single_cell/figure2/number-of-expressed-receptors-cellTypeGroup.pdf', width = 4, height = 2.5)


#################################
# Number of expressed cell types for each HR
df_receptor <- df_receptor[str_order(df_receptor$Gene_symbol),]
df_count <- df_celltype %>% group_by(Gene.name) %>% summarize(count = sum(nTPM > 0)) %>% arrange(count)
df_count$Gene.name <- factor(df_count$Gene.name, levels = df_count$Gene.name)
df_count <- df_count[str_order(df_count$Gene.name),]
df_count$Receptor_subcellular <- df_receptor$Receptor_subcellular

mycoms <- list(c('Membrane', 'Nucleus'))
ggplot(df_count, aes(Receptor_subcellular, count)) +
  geom_boxplot(aes(fill = Receptor_subcellular), outlier.shape = NA) +
  #geom_jitter(width = 0.3) +
  theme_classic() +
  scale_fill_manual(values = c("#ED1B28","#237EB6"))+
  theme(legend.position = 'none') +
  stat_compare_means(method = 'wilcox.test', comparisons = mycoms) +
  labs(x=NULL, y='Number of cell types\nexpressing a certain HR')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
# ggsave('output/HPA-single_cell/figure2/number-of-expressed-cell-types-subcellular.pdf', width = 1.5, height = 2.5)

#################################
# # plot histogram for cut-off：大于40定义为ubiquitously expressed HRs，小于40定义为selectively expressed HRs。
# ggplot(df_count, aes(count)) +
#   geom_density(kernel = 'gaussian') +
#   theme_classic() +
#   labs(y='Distribution density', x='Number of cell types\nexpressing a certain HR')
# # ggsave('output/HPA-single_cell/figure2/number-of-expressed-cell-types-density1.pdf', width = 2.5, height = 2.5)
# 
# # encode HR based on expressed cell types count
# df_count$Expression_universality <- to_vec(for(i in df_count$count) if(i > 40) 'ubiquitously expressed HRs' else 'selectively expressed HRs')
# #write_tsv(df_count, 'input/expression_universality.tsv')
# ggplot(df_count, aes(Expression_universality, fill=Receptor_subcellular)) +
#   geom_bar(position='stack') +
#   scale_fill_manual(values = c("#ED1B28","#237EB6"))+
#   theme_classic() +
#   labs(x=NULL, y='Count',fill='HR types')+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
# # ggsave('output/HPA-single_cell/figure2/expression_universality.pdf', width = 2.5, height = 3)

#################################
# plot histogram for cut-off：分n类，cutoff分别是总细胞数的10%（8），30%（24），70%（55），90%（71）
ggplot(df_count, aes(count)) +
  geom_density(kernel = 'gaussian') +
  geom_vline(xintercept=8,linetype = "dashed")+
  geom_vline(xintercept=24,linetype = "dashed")+
  geom_vline(xintercept=55,linetype = "dashed")+
  geom_vline(xintercept=71,linetype = "dashed")+
  theme_classic() +
  labs(y='Distribution density', x='Number of cell types\nexpressing a certain HR')
# ggsave('output/HPA-single_cell/figure2/number-of-expressed-cell-types-density-new.pdf', width = 2.5, height = 2.5)

# encode HR based on expressed cell types count
df_count$Expression_universality <- to_vec(for(i in df_count$count) if(i<8) '< 10%' else if(i<24) '10% - 30%' else if(i<55) '30% - 70%' else if(i<71) '70% - 90%' else '> 90%')
#write_tsv(df_count, 'input/expression_universality.tsv')
df_count$Expression_universality <- factor(df_count$Expression_universality, levels = c('< 10%','10% - 30%','30% - 70%','70% - 90%','> 90%'))
ggplot(df_count, aes(Expression_universality, fill=Receptor_subcellular)) +
  geom_bar(position='stack') +
  scale_fill_manual(values = c("#ED1B28","#237EB6"))+
  theme_classic() +
  labs(x='Percentage of total cell type number', y='Count',fill='HR types')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
# ggsave('output/HPA-single_cell/figure2/expression_universality-new.pdf', width = 3.5, height = 2.5)


# encode HR based on expressed cell types count -- 1
df_count$Expression_universality <- to_vec(for(i in df_count$count) if(i<10) '< 10' else if(i<20) '10 - 20' else if(i<40) '20 - 40' else if(i<60) '40 - 60' else if(i<70) '60 - 70' else '> 70')
#write_tsv(df_count, 'input/expression_universality.tsv')
df_count$Expression_universality <- factor(df_count$Expression_universality, levels = c('< 10','10 - 20','20 - 40','40 - 60','60 - 70','> 70'))
ggplot(df_count, aes(Expression_universality, fill=Receptor_subcellular)) +
  geom_bar(position='stack') +
  scale_fill_manual(values = c("#ED1B28","#237EB6"))+
  theme_classic() +
  labs(x='No. of expressed cell types', y='No. of HRs',fill='HR types')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black")) 
# ggsave('output/HPA-single_cell/figure2/expression_universality-new-1.pdf', width = 3, height = 2.5)


# --------------------------
# print the specific gene names in each range.
g <- df_count %>% filter(Expression_universality == '< 10%') %>% pull(Gene.name) # < 10%
# CCKAR CRHR1 DRD5  MC2R  MC3R  MC5R  MLNR  NPY4R RXFP2 RXFP3 SSTR4 SSTR5
df <- df_celltype %>% filter(Gene.name %in% g)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
cell.type.tokeep <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% filter(n > 0) %>% pull(Cell.type)
df <- df %>% filter(Cell.type %in% cell.type.tokeep)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name))+
  geom_bar(stat = 'identity', position = 'stack')+
  coord_flip()+
  scale_fill_manual(values = col_vector) +
  labs(x='Cell type', fill='Gene name') +
  theme_classic()
# ggsave('output/HPA-single_cell/figure2/expression_universality-new-10percentage.pdf', width = 5, height = 5)

g <- df_count %>% filter(Expression_universality == '10% - 30%') %>% pull(Gene.name) # 10% - 30%
# ADRA1D ADRB3  AGTR2  AVPR1B CALCR  CRHR2  DRD1   DRD3   FSHR   GHRHR  GHSR   GNRHR  HCRTR1 HCRTR2 LHCGR  MC4R MTNR1A NR2E1  NR5A1  RXFP4  SSTR3  TRHR
df <- df_celltype %>% filter(Gene.name %in% g)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
cell.type.tokeep <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% filter(n > 0) %>% pull(Cell.type)
df <- df %>% filter(Cell.type %in% cell.type.tokeep)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name))+
  geom_bar(stat = 'identity', position = 'stack')+
  coord_flip()+
  scale_fill_manual(values = col_vector) +
  labs(x='Cell type', fill='Gene name') +
  theme_classic()
# ggsave('output/HPA-single_cell/figure2/expression_universality-new-10percentageTo30percentage.pdf', width = 6, height = 10)


g <- df_count %>% filter(Expression_universality == '70% - 90%') %>% pull(Gene.name) # 70% - 90%
df <- df_celltype %>% filter(Gene.name %in% g)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
cell.type.tokeep <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% filter(n > 0) %>% pull(Cell.type)
df <- df %>% filter(Cell.type %in% cell.type.tokeep)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name))+
  geom_bar(stat = 'identity', position = 'stack')+
  coord_flip()+
  scale_fill_manual(values = col_vector) +
  labs(x='Cell type', fill='Gene name') +
  theme_classic()
# ggsave('output/HPA-single_cell/figure2/expression_universality-new-70percentageTo90percentage.pdf', width = 8, height = 10)

g <- df_count %>% filter(Expression_universality == '> 90%') %>% pull(Gene.name) # > 90%
df <- df_celltype %>% filter(Gene.name %in% g)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
cell.type.tokeep <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% filter(n > 0) %>% pull(Cell.type)
df <- df %>% filter(Cell.type %in% cell.type.tokeep)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name))+
  geom_bar(stat = 'identity', position = 'stack')+
  coord_flip()+
  scale_fill_manual(values = col_vector) +
  labs(x='Cell type', fill='Gene name') +
  theme_classic()
# ggsave('output/HPA-single_cell/figure2/expression_universality-new-90percentage.pdf', width = 9, height = 10)

#################################
# HR module
df_module <- readxl::read_excel('input/HR-module.xlsx')
df_module <- left_join(df_module, df_receptor, 'Gene_symbol') %>% filter(Hormone_name1 != 'orphan')
ggplot(df_module, aes(Module, fill=Hormone_name1)) +
  geom_bar(position = 'fill') +
  coord_flip() +
  scale_fill_manual(values = col_vector)

df <- df_module %>% group_by(Hormone_name1,Hormone_chemical_classes) %>% count()
set.seed(42)
ggplot(df, aes(label = Hormone_name1, size = n, color=Hormone_chemical_classes)) +
  geom_text_wordcloud() +
  theme_minimal() +
  scale_color_manual(values = c('#87A0C7','#FF8D6C','#41C2A6'))
# ggsave('output/HPA-single_cell/figure2/HR-module-wordCloud.pdf', width = 2, height = 2)


df_123 <- df_module %>% filter(Module %in% c('Module 1','Module 2','Module 3')) %>% group_by(Hormone_name1,Hormone_chemical_classes) %>% count()
set.seed(42)
ggplot(df_123, aes(label = Hormone_name1, size = n, color=Hormone_chemical_classes)) +
  geom_text_wordcloud() +
  theme_minimal() +
  scale_color_manual(values = c('#87A0C7','#FF8D6C','#41C2A6'))
# ggsave('output/HPA-single_cell/figure2/HR-module-wordCloud-module123.pdf', width = 2, height = 2)

df_4 <- df_module %>% filter(Module == 'Module 4') %>% group_by(Hormone_name1,Hormone_chemical_classes) %>% count()
set.seed(42)
ggplot(df_4, aes(label = Hormone_name1, size = n, color=Hormone_chemical_classes)) +
  geom_text_wordcloud() +
  theme_minimal() +
  scale_color_manual(values = c('#87A0C7','#FF8D6C','#41C2A6'))
# ggsave('output/HPA-single_cell/figure2/HR-module-wordCloud-module4.pdf', width = 2, height = 2)


#################################
# Expression analysis of adrenergic receptors
vec_ADRs <- c('ADRA1A','ADRA1B','ADRA1D','ADRA2A','ADRA2B','ADRA2C','ADRB1','ADRB2','ADRB3')
df <- df_celltype %>% filter(Gene.name %in% vec_ADRs)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name)) +
  geom_bar(stat='identity', position = 'stack')+
  coord_flip()+
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  labs(x='Cell type', fill='Gene name')
# ggsave('output/HPA-single_cell/figure2/HR-module-ADR.pdf', width = 5, height = 9)

# Expression analysis of dopamine receptors
vec_DRD <- c('DRD1','DRD2','DRD3','DRD4','DRD5')
df <- df_celltype %>% filter(Gene.name %in% vec_DRD)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name)) +
  geom_bar(stat='identity', position = 'stack')+
  coord_flip()+
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  labs(x='Cell type', fill='Gene name')
# ggsave('output/HPA-single_cell/figure2/HR-module-DRD.pdf', width = 5, height = 9)

# Expression analysis of proglucagon (GCG)-related receptors
vec_GCG <- c('GCGR','GLP1R','GLP2R')
df <- df_celltype %>% filter(Gene.name %in% vec_GCG)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name)) +
  geom_bar(stat='identity', position = 'stack')+
  coord_flip()+
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  labs(x='Cell type', fill='Gene name')
# ggsave('output/HPA-single_cell/figure2/HR-module-GCG.pdf', width = 5, height = 9)

# Expression analysis of somatostatin receptors
vec_SSTR <- c('SSTR1','SSTR2','SSTR3','SSTR4','SSTR5')
df <- df_celltype %>% filter(Gene.name %in% vec_SSTR)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name)) +
  geom_bar(stat='identity', position = 'stack')+
  coord_flip()+
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  labs(x='Cell type', fill='Gene name')
# ggsave('output/HPA-single_cell/figure2/HR-module-SSTR.pdf', width = 5, height = 9)


# Expression analysis of Nuclear Receptor Subfamily 1
vec_NR1 <- c('NR1D1',
             'NR1D2',
             'NR1H2',
             'NR1H3',
             'NR1H4',
             'NR1I2',
             'NR1I3')
df <- df_celltype %>% filter(Gene.name %in% vec_NR1)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name)) +
  geom_bar(stat='identity', position = 'stack')+
  coord_flip()+
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  labs(x='Cell type', fill='Gene name') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
ggsave('output/HPA-single_cell/figure2/HR-module-NR1.pdf', width = 5, height = 9)

# Expression analysis of Nuclear Receptor Subfamily 2
vec_NR2 <- c('NR2C1',
             'NR2C2',
             'NR2E1',
             'NR2E3',
             'NR2F1',
             'NR2F2',
             'NR2F6')
df <- df_celltype %>% filter(Gene.name %in% vec_NR2)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name)) +
  geom_bar(stat='identity', position = 'stack')+
  coord_flip()+
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  labs(x='Cell type', fill='Gene name') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
ggsave('output/HPA-single_cell/figure2/HR-module-NR2.pdf', width = 5, height = 9)


# Expression analysis of various GPCR families: class A, classB, class C...


# Expression analysis of various nuclear hormone receptor (NHR) families:


#################################
#-----------------top 10-----------------
# Expression analysis of adrenergic receptors
vec_ADRs <- c('ADRA1A','ADRA1B','ADRA1D','ADRA2A','ADRA2B','ADRA2C','ADRB1','ADRB2','ADRB3')
df <- df_celltype %>% filter(Gene.name %in% vec_ADRs)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
df <- df %>% filter(Cell.type %in% rev(celltype.order)[1:10])
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name)) +
  geom_bar(stat='identity', position = 'stack')+
  coord_flip()+
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  labs(x='Cell type', fill='Gene name') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure2/HR-module-ADR-top10.pdf', width = 5, height = 2.5)

# Expression analysis of dopamine receptors
vec_DRD <- c('DRD1','DRD2','DRD3','DRD4','DRD5')
df <- df_celltype %>% filter(Gene.name %in% vec_DRD)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
df <- df %>% filter(Cell.type %in% rev(celltype.order)[1:10])
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name)) +
  geom_bar(stat='identity', position = 'stack')+
  coord_flip()+
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  labs(x='Cell type', fill='Gene name') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure2/HR-module-DRD-top10.pdf', width = 4.6, height = 2.5)

# Expression analysis of proglucagon (GCG)-related receptors
vec_GCG <- c('GCGR','GLP1R','GLP2R')
df <- df_celltype %>% filter(Gene.name %in% vec_GCG)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
df <- df %>% filter(Cell.type %in% rev(celltype.order)[1:10])
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name)) +
  geom_bar(stat='identity', position = 'stack')+
  coord_flip()+
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  labs(x='Cell type', fill='Gene name') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure2/HR-module-GCG-top10.pdf', width = 4.6, height = 2.5)

# Expression analysis of somatostatin receptors
vec_SSTR <- c('SSTR1','SSTR2','SSTR3','SSTR4','SSTR5')
df <- df_celltype %>% filter(Gene.name %in% vec_SSTR)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
df <- df %>% filter(Cell.type %in% rev(celltype.order)[1:10])
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name)) +
  geom_bar(stat='identity', position = 'stack')+
  coord_flip()+
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  labs(x='Cell type', fill='Gene name') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure2/HR-module-SSTR-top10.pdf', width = 5, height = 2.5)

# Expression analysis of Nuclear Receptor Subfamily 1
vec_NR1 <- c('NR1D1',
             'NR1D2',
             'NR1H2',
             'NR1H3',
             'NR1H4',
             'NR1I2',
             'NR1I3')
df <- df_celltype %>% filter(Gene.name %in% vec_NR1)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
df <- df %>% filter(Cell.type %in% rev(celltype.order)[1:10])
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name)) +
  geom_bar(stat='identity', position = 'stack')+
  coord_flip()+
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  labs(x='Cell type', fill='Gene name') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure2/HR-module-NR1-top10.pdf', width = 4.7, height = 2.5)

# Expression analysis of Nuclear Receptor Subfamily 2
vec_NR2 <- c('NR2C1',
             'NR2C2',
             'NR2E1',
             'NR2E3',
             'NR2F1',
             'NR2F2',
             'NR2F6')
df <- df_celltype %>% filter(Gene.name %in% vec_NR2)
celltype.order <- df %>% group_by(Cell.type) %>% summarise(n=sum(nTPM)) %>% arrange(n) %>% pull(Cell.type)
df$Cell.type <- factor(df$Cell.type, levels = celltype.order)
df <- df %>% filter(Cell.type %in% rev(celltype.order)[1:10])
ggplot(df, aes(Cell.type, nTPM, fill=Gene.name)) +
  geom_bar(stat='identity', position = 'stack')+
  coord_flip()+
  theme_classic() +
  scale_fill_brewer(palette = 'Set1') +
  labs(x='Cell type', fill='Gene name') +
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
# ggsave('output/HPA-single_cell/figure2/HR-module-NR2-top10.pdf', width = 4.6, height = 2.5)


# Expression analysis of various GPCR families: class A, classB, class C...


# Expression analysis of various nuclear hormone receptor (NHR) families:


