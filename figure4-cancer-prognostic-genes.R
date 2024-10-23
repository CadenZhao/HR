#################################
# Author: Xiangjie Zhao
# Date: 2023.04.20
# Title: Analysis of cancer prognostic genes
#################################

rm(list = ls())
set.seed(1)
library(tidyverse)
library(ggpie)
library(ggstatsplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggvenn)
setwd('~/project/hormone_receptor_scRNA-2022-12-11/')

#################################
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
df_receptor$Hormone_chemical_classes <- recode(df_receptor$Hormone_chemical_classes, 'NA'='Orphan')
vec_receptor <- df_receptor$Gene_symbol

df_cancer <- readxl::read_excel('input/pathology.xlsx')
df_cancer <- df_cancer %>% filter(`Gene name` %in% vec_receptor)
df_cancer <- df_cancer %>% pivot_longer(cols = c("prognostic - favorable","unprognostic - favorable","prognostic - unfavorable","unprognostic - unfavorable"),
                                        names_to = 'Prognosis', values_to = 'p-value')
df_cancer <- df_cancer %>% drop_na('p-value')
df_cancer <- df_cancer %>% separate(Prognosis, sep = " - ", into = c("prognosis", "clinical_outcome"))
# df_cancer$Cancer <- factor(df_cancer$Cancer, levels = rev(df_cancer$Cancer %>% unique()))
df_cancer$Cancer <- factor(df_cancer$Cancer, levels = df_cancer$Cancer %>% unique())

df_cancer_prognostic <- df_cancer %>% filter(prognosis == 'prognostic')
df_cancer_unprognostic <- df_cancer %>% filter(prognosis == 'unprognostic')
# df_cancer_prognostic %>% write_tsv('output/HPA-single_cell/figure4/prognosticHR-clinicalOutcome.tsv')

#################################
# pie plot: Prognostic and Unprognostic
receptor_prognostic <- df_cancer_prognostic$`Gene name` %>% unique()
n_receptor_prognostic <- receptor_prognostic %>% length()
n_receptor_unprognostic <- length(vec_receptor) - n_receptor_prognostic
label <- c('Prognostic','Unprognostic')
n <- c(n_receptor_prognostic,n_receptor_unprognostic)
percentage <- round(100 * n / sum(n), 0)
labs <- paste0(label, " ", " (", percentage, "%)")
# pdf('output/HPA-single_cell/figure4/pie_prognosis.pdf', width = 4, height = 4)
pie(c(n_receptor_prognostic,n_receptor_unprognostic), labels = labs, col=c('#57C2A6','#999999'))
# dev.off()

#################################
# write prognostic and unprognostic HRs.
receptor_unprognostic <- setdiff(vec_receptor, receptor_prognostic)
df1 <- data.frame(gene = receptor_prognostic, prognosis = 'prognostic')
df2 <- data.frame(gene = receptor_unprognostic, prognosis = 'unprognostic')
df12 <- rbind(df1, df2)
# df12 %>% write_tsv('output/HPA-single_cell/figure4/HR_prognosis.tsv')

#################################
df12 <- as_tibble(df12) %>% dplyr::rename(Gene_symbol = gene)
df <- left_join(df_receptor, df12, by='Gene_symbol')

ggbarstats(df, prognosis, Receptor_subcellular, results.subtitle=F)+
  theme_classic()+
  labs(x='HR types', y='Percentage')+
  scale_fill_manual(values = c("#999999","#74C0A7"))
# ggsave('output/HPA-single_cell/figure4/prognosis-HR_types.pdf', width = 3.5, height = 3)

ggbarstats(df, prognosis, Hormone_chemical_classes, results.subtitle=F)+
  theme_classic()+
  labs(x='Hormone classes', y='Percentage')+
  scale_fill_manual(values = c("#999999","#74C0A7"))
# ggsave('output/HPA-single_cell/figure4/prognosis-Hormone_class.pdf', width = 5, height = 3)


df$prognosisReceptor_subcellular <- str_c(df$prognosis, '-', df$Receptor_subcellular)
ggbarstats(df, prognosisReceptor_subcellular, Hormone_chemical_classes, results.subtitle=F)+
  theme_classic()+
  labs(x='Hormone classes', y='Percentage')+
  scale_fill_brewer(palette = 'Paired')
# ggsave('output/HPA-single_cell/figure4/prognosis-Hormone_class-new.pdf', width = 6, height = 3)


#################################
# bar plot: favorable and unfavorable
df_cancer_prognostic$clinical_outcome <- factor(df_cancer_prognostic$clinical_outcome, levels = c('unfavorable','favorable'))
ggplot(df_cancer_prognostic, aes(Cancer, fill=clinical_outcome))+
  geom_bar(position = 'dodge', width = 0.6) +
  labs(y='Count', x='Cancer types', fill='Clinical outcome') +
  # coord_flip() +
  theme_classic() +
  scale_fill_manual(values = c("#92010A","#008C1A")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(color="black"), axis.ticks = element_line(color = "black"))
ggsave('output/HPA-single_cell/figure4/clinical_outcome-1.pdf', width = 5, height = 2.5)


# venn diagram: favorable in all cancers, unfavorable in all cancers
# and the intersection representing HRs can be favorable in some cancer but be unfavorable in other cancers
HR_favorable <- df_cancer_prognostic %>% filter(clinical_outcome == 'favorable') %>% pull('Gene name')
HR_unfavorable <- df_cancer_prognostic %>% filter(clinical_outcome == 'unfavorable') %>% pull('Gene name')
x <- list( `Unfavorable in all cancers` = HR_unfavorable, `Favorable in all cancers` = HR_favorable)
ggvenn(x, fill_color = c("#92010A","#008C1A"), fill_alpha=0.7,
  stroke_size = 0.5, set_name_size = 4)
# ggsave('output/HPA-single_cell/figure4/venn-clinical_outcome.pdf', width = 3.5, height = 2.5)

HR_allFavorable <- setdiff(HR_favorable, intersect(HR_favorable,HR_unfavorable))
HR_allUnfavorable <- setdiff(HR_unfavorable, intersect(HR_favorable,HR_unfavorable))

df1 <- data.frame(gene = HR_allFavorable, clinical_outcome = 'allFavorable')
df2 <- data.frame(gene = HR_allUnfavorable, clinical_outcome = 'allUnfavorable')
df3 <- data.frame(gene = intersect(HR_favorable,HR_unfavorable), clinical_outcome = 'both')
df123 <- rbind(df1, df2, df3)
df123 %>% write_tsv('output/HPA-single_cell/figure4/allAndBoth-clinicalOutcome.tsv')
