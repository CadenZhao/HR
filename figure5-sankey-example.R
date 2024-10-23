#################################
# Author: Xiangjie Zhao
# Date: 2023.05.25
# Title: the organ-organ communication of single hormone level
#################################

rm(list = ls())
set.seed(1)
library(tidyverse)
library(ggsankey)
library(RColorBrewer)
library(comprehenr)
setwd('~/project/hormone_receptor_scRNA-2022-12-11/')

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
male_tissue <- c('Testis', 'Prostate')
female_tissue <- c('Fallopian Tube', 'Vagina', 'Ovary', 'Uterus', 'Cervix Uteri')
node_order_male <- c("Peptide","Lipid-derived","Amino acid-derived","Nucleus HR","Membrane HR","Testis","Prostate","Thyroid","Stomach","Spleen","Small Intestine","Skin","Salivary Gland","Pituitary","Pancreas","Nerve","Muscle","Lung","Liver","Kidney","Heart","Esophagus","Colon","Breast","Brain","Blood Vessel","Blood","Bladder","Adrenal Gland","Adipose Tissue")
node_order_male_color <- c(col_vector[2:3],col_vector[30],col_vector[5:6],col_vector[68:69],col_vector[7:29])
names(node_order_male_color) <- node_order_male
node_order_female <- c("Peptide","Lipid-derived","Amino acid-derived","Nucleus HR","Membrane HR","Vagina","Uterus","Ovary","Fallopian Tube","Cervix Uteri","Thyroid","Stomach","Spleen","Small Intestine","Skin","Salivary Gland","Pituitary","Pancreas","Nerve","Muscle","Lung","Liver","Kidney","Heart","Esophagus","Colon","Breast","Brain","Blood Vessel","Blood","Bladder","Adrenal Gland","Adipose Tissue")
node_order_female_color <- c(col_vector[2:3],col_vector[30],col_vector[5:6],col_vector[70:74],col_vector[7:29])
names(node_order_female_color) <- node_order_female

##################################################################
# 读取激素和HR的信息
##################################################################
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
df_receptor$Hormone_name1 <- recode(df_receptor$Hormone_name1, 'orphan'='Unknown')
df_receptor$Hormone_chemical_classes <- recode(df_receptor$Hormone_chemical_classes, 'NA'='Unknown')
df_receptor$Hormone1_tissue <- recode(df_receptor$Hormone1_tissue, 'NA'='Unknown')
vec_receptor <- df_receptor$Gene_symbol
df_receptor <- df_receptor %>% separate_rows(Hormone1_tissue, sep = ',')
df_receptor <- df_receptor %>% filter(Hormone1_tissue != 'Unknown') # 分析组织通讯关系时候不考虑孤儿受体

df_receptor$Receptor_subcellular <- df_receptor$Receptor_subcellular %>% recode(Membrane = 'Membrane HR', Nucleus = 'Nucleus HR')


##################################################################
# 计算HR表达的细胞和组织：选表达前2高的细胞/组织 --- 分性别
##################################################################
# --------------男性--------------
# 计算最高表达每个HR的细胞类型（HPA数据，79种）
df_celltype <- read.table('input/HPA-single_cell/single-cell-type/rna_single_cell_type.tsv',
                          sep = '\t', header = T)
df_celltype <- filter(df_celltype, Gene.name %in% vec_receptor)
df_celltype$Cell.type <- df_celltype$Cell.type %>% str_to_title
df_celltype <- df_celltype %>% rename('Gene_symbol'='Gene.name')
df_celltype <- df_celltype %>% group_by(Gene_symbol) %>% filter(nTPM > 0, nTPM >= sort(nTPM, TRUE)[2])
df <- left_join(df_receptor, df_celltype, by='Gene_symbol')
df <- df %>% filter(!Hormone1_tissue %in% female_tissue) # 只留该性别有的分泌激素的组织
# 计算最高表达每个HR的组织类型（GTEx数据，30种）
gtex_hr <- read_tsv('input/gtex/GTEx_HR_combined.tsv')
gtex_hr <- filter(gtex_hr, Description %in% vec_receptor, SEX == 'Male') # 筛选性别
gtex_hr <- gtex_hr %>% group_by(Description, SMTS) %>% summarise(TPM=mean(TPM))
gtex_hr <- gtex_hr %>% rename('Gene_symbol'='Description')
gtex_hr <- gtex_hr %>% group_by(Gene_symbol) %>% filter(log1p(TPM) > 0.1, TPM >= sort(TPM, TRUE)[2]) # 此阈值的选定见Fig. S2B
df_merged <- left_join(df, gtex_hr, by='Gene_symbol')
df_merged <- df_merged %>% add_row(Hormone1_tissue = c("Spleen","Esophagus","Prostate","Muscle","Breast" ,"Bladder")) # 补全激素组织
df_merged <- df_merged %>% rename('Hormone name' = 'Hormone_name1',
                                  'Hormone class' = 'Hormone_chemical_classes',
                                  'Hormone tissue' = 'Hormone1_tissue',
                                  'HR tissue' = 'SMTS',
                                  'HR type' = 'Receptor_subcellular',
                                  'HR cell type' = 'Cell.type',
                                  'HR gene name' = 'Gene_symbol')
# sankey plot：summary（组织在中间）
df_merged <- df_merged %>% select(-`HR cell type`) # 只统计器官通讯时，去掉HR细胞的信息

###########################
# glucagon
df <- df_merged %>% filter(`Hormone name` == 'glucagon') %>% make_long(`Hormone class`, `Hormone tissue`, `HR tissue`, `HR type`)
df$node <- factor(df$node, levels = node_order_male) # 按照性别通用组织排序，然后再排性别特异的组织
df$next_node <- factor(df$next_node, levels = node_order_male) # 按照性别通用组织排序，然后再排性别特异的组织
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
  geom_sankey() +
  geom_sankey_label() +
  theme_sankey() +
  theme(legend.position = "none")+
  scale_fill_manual(values = node_order_male_color) +
  labs(x = NULL)
ggsave('output/HPA-single_cell/figure5/sankey-example-glucagon.pdf', width = 4.5, height = 3)

###########################
# CMKLR1
df <- df_merged %>% filter(`HR gene name` == 'CMKLR1') %>% make_long(`Hormone class`, `Hormone tissue`, `HR tissue`, `HR type`)
df$node <- factor(df$node, levels = node_order_male) # 按照性别通用组织排序，然后再排性别特异的组织
df$next_node <- factor(df$next_node, levels = node_order_male) # 按照性别通用组织排序，然后再排性别特异的组织
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
  geom_sankey() +
  geom_sankey_label() +
  theme_sankey() +
  theme(legend.position = "none")+
  scale_fill_manual(values = node_order_male_color) +
  labs(x = NULL)
# ggsave('output/HPA-single_cell/figure5/sankey-example-CMKLR1.pdf', width = 5, height = 2.5)

###########################
# HR tissue: brain
df <- df_merged %>% filter(`HR tissue` == 'Brain') %>% make_long(`Hormone class`, `Hormone tissue`, `HR tissue`, `HR type`)
df$node <- factor(df$node, levels = node_order_male) # 按照性别通用组织排序，然后再排性别特异的组织
df$next_node <- factor(df$next_node, levels = node_order_male) # 按照性别通用组织排序，然后再排性别特异的组织
ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
  geom_sankey() +
  geom_sankey_label() +
  theme_sankey() +
  theme(legend.position = "none")+
  scale_fill_manual(values = node_order_male_color) +
  labs(x = NULL)
ggsave('output/HPA-single_cell/figure5/sankey-example-2brain.pdf', width = 5.7, height = 4)





