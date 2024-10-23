#################################
# Author: Xiangjie Zhao
# Date: 2022.12.24
# Title: Analysis of disease-gene data
#################################

rm(list = ls())
set.seed(1)
library(tidyverse)
library(ggpie)
setwd('~/project/hormone_receptor_scRNA-2022-12-11/')

#################################
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
df_receptor$Hormone_chemical_classes <- recode(df_receptor$Hormone_chemical_classes, 'NA'='Orphan')

ggpie::ggpie(data = df_receptor, group_key = "Receptor_subcellular", count_type = "full",
      label_info = "all", label_type = "circle", label_split = NULL, label_color = "white",
      border_color = "white", label_pos = "in") +
  scale_fill_manual(values = c("#ED1B28","#237EB6"))

# ggsave('output/HPA-single_cell/figure1/pie_subcellular1.pdf', width = 4, height = 4)
  
#################################
# HR subcellular location classified by hormone chemical classes (HPA annotation based on protein IF experiments)
df_subcellular <- readxl::read_excel('input/subcellular_location.xlsx') %>% dplyr::select(`Gene name`, Reliability, `Main location`)
df_subcellular$Gene_symbol <- df_subcellular$`Gene name`
df <- left_join(df_receptor, df_subcellular, by='Gene_symbol')
df <- df %>% mutate(`Main location` = ifelse(is.na(`Main location`), 'Unannotated', `Main location`))
df <- df %>% separate_rows(`Main location`, sep = ';')

df$`Main location` <- factor(df$`Main location`, levels = sort(unique(df$`Main location`), decreasing = T),)


vec_Cytosol <- c('Cytoplasmic bodies','Cytosol','Focal adhesion sites','Golgi apparatus','Microtubules')
vec_Nucleoli <- c('Nucleoli','Nucleoli fibrillar center')
vec_Nucleoplasm <- c('Nucleoplasm','Nuclear speckles','Nuclear bodies')
vec_Cell_Junctions <- c('Cell Junctions')
vec_Nuclear_membrane <- c('Nuclear membrane')
vec_Plasma_membrane <- c('Plasma membrane')
vec_Unannotated <- c('Unannotated')
vec_Vesicles <- c('Vesicles')

df <- df %>% mutate(`Merged main location` = case_when(`Main location` %in% vec_Cytosol ~ 'Cytosol',
                                                       `Main location` %in% vec_Nucleoli ~ 'Nucleoli',
                                                       `Main location` %in% vec_Nucleoplasm ~ 'Nucleoplasm',
                                                       `Main location` %in% vec_Cell_Junctions ~ 'Cell Junctions',
                                                       `Main location` %in% vec_Nuclear_membrane ~ 'Nuclear membrane',
                                                       `Main location` %in% vec_Plasma_membrane ~ 'Plasma membrane',
                                                       `Main location` %in% vec_Unannotated ~ 'Unannotated',
                                                       `Main location` %in% vec_Vesicles ~ 'Vesicles',
                                                       ))

df$`Merged main location` <- factor(df$`Merged main location`, levels = sort(unique(df$`Merged main location`), decreasing = T))

ggplot(df, aes(`Merged main location`, fill=Hormone_chemical_classes))+
  geom_bar() +
  theme_classic()+
  labs(x=NULL, y='Number of genes', fill='Hormone classes')+
  # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  coord_flip()+
  theme(axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
ggsave('output/HPA-single_cell/figure1/HR-subcellularLocation-by-hormone-class-Merged.pdf', width = 4.5, height = 4)



