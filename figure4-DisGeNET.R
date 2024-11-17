#################################
# Author: Xiangjie Zhao
# Date: 2022.12.24
# Title: Analysis of disease-gene data
#################################

rm(list = ls())
set.seed(1)
library(tidyverse)
library(disgenet2r)
library(ggpie)
setwd('~/project/hormone_receptor_scRNA-2022-12-11/')

#################################
disgenet_api_key <- get_disgenet_api_key(
  email = "**********", 
  password = "**********" )
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

#################################
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
vec_receptor <- df_receptor$Gene_symbol

dis_receptor <- gene2disease(gene = vec_receptor,
                             score = c(0, 1),
                             database = "ALL")
dis_receptor_qresult <- dis_receptor@qresult
dis_receptor_qresult["disease_class_name"][dis_receptor_qresult["disease_class_name"] == ''] <- 'Other'
dis_receptor_qresult["disease_name"][dis_receptor_qresult["disease_name"] == ''] <- 'Other'
dis_receptor_qresult["protein_class_name"][dis_receptor_qresult["protein_class_name"] == ''] <- 'Other'

#################################
# disease class name
# select top associated diseases for each receptor
dis_receptor_top <- (dis_receptor_qresult %>% as_tibble %>% group_by(gene_symbol) %>% top_n(1, score)
                     %>% arrange(gene_symbol) %>% group_by(gene_symbol) %>% dplyr::slice(1)
                     %>% separate(disease_class_name, c('disease_class_name_1','disease_class_name_2'), sep = ';')
                     )
dis_receptor_top$disease_class_name_1 <- str_trim(dis_receptor_top$disease_class_name_1)

#################################
# save and manually refine the disease class

# write.table(dis_receptor_top, file = 'output/HPA-single_cell/table1/disease-pie1.tsv',
#              quote = F, sep = '\t', row.names = F)
dis_receptor_top <- readxl::read_excel('output/HPA-single_cell/table1/disease-pie1.xlsx')

#################################
# plot
ggrosepie(dis_receptor_top, c("disease_class_name_1","protein_class_name"),
          count_type = "full", label_info = "all", label_sep = ' #',
          donut_frac=0.3, donut_label = F, label_size = 7)
# ggsave('output/HPA-single_cell/figure1/disease-pie1.pdf', width = 8, height = 6)

# plot disease name horizontally
ggplot(dis_receptor_top, aes(protein_class_name, fill=disease_class_name_1)) +
  geom_bar()
# ggsave('output/HPA-single_cell/figure1/disease-pie1-name.pdf', width = 6, height = 6)

