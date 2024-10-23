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
library(RColorBrewer)

#################################
# generate a number of most distinctive colors in R?
# 74 colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


setwd('~/project/hormone_receptor_scRNA-2022-12-11/')
color=rev(brewer.pal(11, name = 'RdYlBu')) # define a heatmap color

#################################
# HR annotations
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
df_receptor <- df_receptor %>% arrange(Receptor_subcellular, Gene_symbol) # fix a gene order to present
vec_receptor <- df_receptor$Gene_symbol

#################################
# Gene - cell.type expression data
df_celltype <- read.table('input/HPA-single_cell/single-cell-type/rna_single_cell_type.tsv',
                          sep = '\t', header = T)
df_celltype <- filter(df_celltype, Gene.name %in% vec_receptor)
df_celltype$Cell.type <- df_celltype$Cell.type %>% str_to_title

#################################
# Format expression data to a matrix
df_celltype_wider <- df_celltype %>% pivot_wider(names_from = Cell.type, values_from = nTPM)
row_name <- df_celltype_wider$Gene.name
df_celltype_wider <- df_celltype_wider %>% select(-Gene, -Gene.name) %>% as.matrix
rownames(df_celltype_wider) <- row_name
df_celltype_wider_log1p <- t(log1p(df_celltype_wider))
matrix_CMKLR2 <- matrix(0, nrow = dim(df_celltype_wider_log1p)[1]) # 注意，没有CMKLR2基因的表达，手动填0
colnames(matrix_CMKLR2) <- 'CMKLR2'
df_celltype_wider_log1p <- cbind(df_celltype_wider_log1p, matrix_CMKLR2)
df_celltype_wider_log1p <- df_celltype_wider_log1p[,vec_receptor]

vec_celltype <- rownames(df_celltype_wider_log1p) # 79 cell types
#################################
# complexheatmap
#################################
# Gene annotation
# receptor_subcellular
colannot <- df_receptor %>% select(Gene_symbol,Receptor_subcellular) %>% as.data.frame
rownames(colannot) <- colannot$Gene_symbol
colannot <- subset(colannot, select = -Gene_symbol)
column_ha <- HeatmapAnnotation(df = colannot, show_legend=T, show_annotation_name = T)

ht_list = Heatmap(df_celltype_wider_log1p, cluster_rows = T, cluster_columns = T, col = color,
                  row_names_side = 'left',
                  top_annotation = column_ha)
ht_list

#################################
# Cell type annotation
# Expression distribution of hormone receptors in each cell types
df_celltype_tissue <- read.table('input/HPA-single_cell/single-cell-type-tissue-cluster/rna_single_cell_type_tissue.tsv',
                                 sep = '\t', header = T)
df_celltype_tissue <- df_celltype_tissue %>% select(Tissue, Cluster, Cell.type) %>% unique
df_celltype_tissue$Tissue <- df_celltype_tissue$Tissue %>% str_to_title
df_celltype_tissue$Cell.type <- df_celltype_tissue$Cell.type %>% str_to_title
# which and how much proportion of a cell type belonging to a tissue?
rowannot <- df_celltype_tissue %>% group_by(Cell.type, Tissue) %>% summarise(n = n()) %>% mutate(freq = n / sum(n)) %>% select(-n)
rowannot <- rowannot %>% pivot_wider(names_from = Tissue, values_from = freq) # convert into a matrix
rowannot[is.na(rowannot)] <- 0
rowannot <- as.data.frame(rowannot)
rownames(rowannot) <- rowannot$Cell.type
rowannot <- subset(rowannot, select = -Cell.type)
rowannot <- rowannot[vec_celltype,]
# library(circlize)
# state_col = c("#FF0000", "#008000", "#C2E105", "#8A91D0", "#CD5C5C", "#808080", "#000000")
# 
anno_width = unit(3, "cm")
row_ha <- rowAnnotation("pct_st" = anno_barplot(rowannot, bar_width = 1, gp = gpar(fill = unique(col_vector)), width = anno_width), show_annotation_name = FALSE)
ht_list = ht_list + row_ha
ht_list

#################################
# Cell type annotation 2 -- cell type description

df_cell_cluster_description <- read.table('input/HPA-single_cell/single-cell-type-tissue-cluster/rna_single_cell_cluster_description.tsv', sep = '\t', header = T)
df_cell_cluster_description$Cell.type <- df_cell_cluster_description$Cell.type %>% str_to_title

df_cell_cluster_description_group <- df_cell_cluster_description %>% select(Cell.type, Cell.type.group) %>% unique() %>% as.data.frame
rownames(df_cell_cluster_description_group) <- df_cell_cluster_description_group$Cell.type
df_cell_cluster_description_group <- subset(df_cell_cluster_description_group, select = -Cell.type)
df_cell_cluster_description_group <- df_cell_cluster_description_group[vec_celltype,] %>% as.data.frame
rownames(df_cell_cluster_description_group) <- vec_celltype
colnames(df_cell_cluster_description_group) <- 'cell.group'
cell.group.unique <- df_cell_cluster_description_group$cell.group %>% unique()
cell.group.color <- setNames(col_vector[1:15], cell.group.unique)
row_ha1 <- rowAnnotation(df = df_cell_cluster_description_group, show_annotation_name = FALSE, show_legend=F, col=list(cell.group=cell.group.color))
ht_list = ht_list + row_ha1
ht_list

#################################
# Cell type annotation 3 -- cell count
# df_cell_cluster_description_cellcount <- df_cell_cluster_description %>% select(Cell.type, Cell.count)





#################################
# add legend, draw and save plot
lgd_list = list(Legend(labels = colnames(rowannot), title = "Tissue", legend_gp = gpar(fill = unique(col_vector)),  direction = "horizontal"),
                Legend(labels = cell.group.unique, title = "Cell type group", legend_gp = gpar(fill = unique(col_vector)),  direction = "horizontal")) # define legend

pdf('output/HPA-single_cell/figure2/overview-1.pdf', width = 30, height = 16)
draw(ht_list, padding = unit(c(2, 2, 20, 2), "mm"), 
     heatmap_legend_list = lgd_list, heatmap_legend_side = "right")
dev.off()


