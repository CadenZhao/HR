#################################
# Author: Xiangjie Zhao
# Date: 2022.12.27
# Title: Gene ontology enrichment
#################################
rm(list = ls())
set.seed(1)

library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(tidyverse)

setwd('~/project/hormone_receptor_scRNA-2022-12-11/')

#################################
df <- read_tsv('output/HPA-single_cell/figure4/allAndBoth-clinicalOutcome.tsv')
allFavorable <- df %>% filter(clinical_outcome == 'allFavorable') %>% pull(gene)
allUnfavorable <- df %>% filter(clinical_outcome == 'allUnfavorable') %>% pull(gene)
both <- df %>% filter(clinical_outcome == 'both') %>% pull(gene)

#################################
# WikiPathways: allFavorable
# gene symbol to NCBI ID
vec_receptor_ncbi <- bitr(allFavorable, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
ewiki <- enrichWP(vec_receptor_ncbi, organism = "Homo sapiens") 
dim(ewiki)
head(ewiki, 3)

ewiki <- ewiki@result
ewiki$pvalue <- -log10(ewiki$p.adjust)
ewiki <- ewiki %>% filter(p.adjust < 0.05)
ewiki$Description <- factor(ewiki$Description, levels=rev(ewiki$Description))
ewiki <- head(ewiki, 5)

ggplot(ewiki, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='WikiPathways',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure4/pathway-wiki-allFavorable.pdf', width = 7, height = 2)

#################################
# WikiPathways: allUnfavorable
# gene symbol to NCBI ID
vec_receptor_ncbi <- bitr(allUnfavorable, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
ewiki <- enrichWP(vec_receptor_ncbi, organism = "Homo sapiens") 
dim(ewiki)
head(ewiki, 3)

ewiki <- ewiki@result
ewiki$pvalue <- -log10(ewiki$p.adjust)
ewiki <- ewiki %>% filter(p.adjust < 0.05)
ewiki$Description <- factor(ewiki$Description, levels=rev(ewiki$Description))
ewiki <- head(ewiki, 5)

ggplot(ewiki, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='WikiPathways',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure4/pathway-wiki-allUnfavorable.pdf', width = 7, height = 2)

#################################
# WikiPathways: both
# gene symbol to NCBI ID
vec_receptor_ncbi <- bitr(both, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
ewiki <- enrichWP(vec_receptor_ncbi, organism = "Homo sapiens") 
dim(ewiki)
head(ewiki, 3)

ewiki <- ewiki@result
ewiki$pvalue <- -log10(ewiki$p.adjust)
ewiki <- ewiki %>% filter(p.adjust < 0.05)
ewiki$Description <- factor(ewiki$Description, levels=rev(ewiki$Description))
ewiki <- head(ewiki, 5)

ggplot(ewiki, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='WikiPathways',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure4/pathway-wiki-both.pdf', width = 7, height = 2)

