#################################
# Author: Xiangjie Zhao
# Date: 2023.03.30
# Title: Analysis of disease-gene data
#################################

rm(list = ls())
set.seed(1)
library(tidyverse)
library(ggpie)
setwd('~/project/hormone_receptor_scRNA-2022-12-11/')

#################################
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
df_receptor$Hormone_chemical_classes <- recode(df_receptor$Hormone_chemical_classes, 'NA'='Unknown')

#################################
# ggpie::ggpie(data = df_receptor, group_key = "Hormone_chemical_classes", count_type = "full",
#              label_info = "all", label_type = "circle", label_split = NULL, label_color = "white",
#              border_color = "white", label_pos = "in")

df_receptor <- df_receptor %>% group_by(Hormone_chemical_classes) %>% summarise(count=n())
df_receptor$percentage <- round(100 * df_receptor$count / sum(df_receptor$count), 0)
labs <- paste0(df_receptor$Hormone_chemical_classes, " ", df_receptor$count, " (", df_receptor$percentage, "%)")
pdf('output/HPA-single_cell/figure1/hormone-class.pdf', width = 3.5, height = 3.5)
pie(df_receptor$count , labels = labs)
dev.off()
