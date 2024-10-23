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
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')
vec_receptor <- df_receptor$Gene_symbol

#################################
# BP
go <- enrichGO(gene = vec_receptor, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
dim(go)

go_head <- head(go,20)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)

ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='GO biological process',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
# ggsave('output/HPA-single_cell/figure1/pathway-BP.pdf', width = 8, height = 5)

# BP: membrane HRs
vec_receptor_type <- df_receptor %>% filter(Receptor_subcellular == 'Membrane') %>% pull(Gene_symbol)
go <- enrichGO(gene = vec_receptor_type, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
dim(go)
go_head <- head(go,10)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)
ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Membrane HRs (GO biological process)',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
# ggsave('output/HPA-single_cell/figure1/pathway-BP-membraneHR.pdf', width = 8, height = 3)

# BP: nucleus HRs
vec_receptor_type <- df_receptor %>% filter(Receptor_subcellular == 'Nucleus') %>% pull(Gene_symbol)
go <- enrichGO(gene = vec_receptor_type, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
dim(go)
go_head <- head(go,10)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)
ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Nucleus HRs (GO biological process)',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
# ggsave('output/HPA-single_cell/figure1/pathway-BP-nucleusHR.pdf', width = 8, height = 3)


# BP: lipid derived HR
vec_receptor_type <- df_receptor %>% filter(Hormone_chemical_classes == 'Lipid derivative') %>% pull(Gene_symbol)
go <- enrichGO(gene = vec_receptor_type, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
dim(go)
go_head <- head(go,10)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)
ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Lipid derived HRs (GO biological process)',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure1/pathway-BP-Lipid-derivative.pdf', width = 6, height = 3)


# BP: amino acid derived HR
vec_receptor_type <- df_receptor %>% filter(Hormone_chemical_classes == 'Amino acid derivative') %>% pull(Gene_symbol)
go <- enrichGO(gene = vec_receptor_type, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
dim(go)
go_head <- head(go,10)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)
ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Amino acid derived HRs (GO biological process)',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure1/pathway-BP-Amino-acid-derivative.pdf', width = 6, height = 3)


# BP: peptide HR
vec_receptor_type <- df_receptor %>% filter(Hormone_chemical_classes == 'Peptide') %>% pull(Gene_symbol)
go <- enrichGO(gene = vec_receptor_type, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
dim(go)
go_head <- head(go,10)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)
ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Peptide HRs (GO biological process)',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure1/pathway-BP-Peptide.pdf', width = 6, height = 3)


# BP: orphan HR
vec_receptor_type <- df_receptor %>% filter(Hormone_chemical_classes == 'NA') %>% pull(Gene_symbol)
go <- enrichGO(gene = vec_receptor_type, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
dim(go)
go_head <- head(go,10)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)
ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Orphan HRs (GO biological process)',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure1/pathway-BP-orphan-HR.pdf', width = 6, height = 3)

#################################
# MF
go <- enrichGO(gene = vec_receptor, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'MF')
dim(go)
head(go,3)

go_head <- head(go,20)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)

ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='GO molecular function',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))

# ggsave('output/HPA-single_cell/figure1/pathway-MF.pdf', width = 6.5, height = 5)

#################################
# CC
go <- enrichGO(gene = vec_receptor, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'CC')
dim(go)
head(go,3)

go_head <- head(go,20)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)

ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='GO cellular component',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))

# ggsave('output/HPA-single_cell/figure1/pathway-CC.pdf', width = 5.5, height = 5)

#################################
# WikiPathways
# gene symbol to NCBI ID
vec_receptor_ncbi <- bitr(vec_receptor, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
ewiki <- enrichWP(vec_receptor_ncbi, organism = "Homo sapiens") 
dim(ewiki)
head(ewiki, 3)

ewiki_head <- head(ewiki,20)
ewiki_head$Description <- factor(ewiki_head$Description, levels=rev(ewiki_head$Description))
ewiki_head$pvalue <- -log10(ewiki_head$p.adjust)

ggplot(ewiki_head, aes(x=Count, y = Description, fill = pvalue)) +
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

# ggsave('output/HPA-single_cell/figure1/pathway-wiki.pdf', width = 7, height = 5)


#################################
# Reactome Pathway
# gene symbol to NCBI ID
library(ReactomePA)
vec_receptor_ncbi <- bitr(vec_receptor, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
ewiki <- enrichPathway(gene=vec_receptor_ncbi, pvalueCutoff = 0.05, readable=TRUE)
dim(ewiki)
head(ewiki, 3)

ewiki_head <- head(ewiki,10)
ewiki_head$Description <- factor(ewiki_head$Description, levels=rev(ewiki_head$Description))
ewiki_head$pvalue <- -log10(ewiki_head$p.adjust)

ggplot(ewiki_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Reactome Pathway',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))

# ggsave('output/HPA-single_cell/figure1/pathway-reactome.pdf', width = 6, height = 5)

# membrane HR
vec_receptor <- df_receptor %>% filter(Receptor_subcellular == 'Membrane') %>% pull(Gene_symbol)
vec_receptor_ncbi <- bitr(vec_receptor, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
ewiki <- enrichPathway(gene=vec_receptor_ncbi, pvalueCutoff = 0.05, readable=TRUE)
dim(ewiki)
head(ewiki, 3)

ewiki_head <- head(ewiki,10)
ewiki_head$Description <- factor(ewiki_head$Description, levels=rev(ewiki_head$Description))
ewiki_head$pvalue <- -log10(ewiki_head$p.adjust)

ggplot(ewiki_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Membrane HRs (Reactome Pathway)',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure1/pathway-reactome-membrane.pdf', width = 6, height = 3)


# nucleus HR
vec_receptor <- df_receptor %>% filter(Receptor_subcellular == 'Nucleus') %>% pull(Gene_symbol)
vec_receptor_ncbi <- bitr(vec_receptor, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
ewiki <- enrichPathway(gene=vec_receptor_ncbi, pvalueCutoff = 0.05, readable=TRUE)
dim(ewiki)
head(ewiki, 3)

ewiki_head <- head(ewiki,10)
ewiki_head$Description <- factor(ewiki_head$Description, levels=rev(ewiki_head$Description))
ewiki_head$pvalue <- -log10(ewiki_head$p.adjust)

ggplot(ewiki_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Nuclear HRs (Reactome Pathway)',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure1/pathway-reactome-nucleus.pdf', width = 6, height = 3)



# Lipid derivative
vec_receptor <- df_receptor %>% filter(Hormone_chemical_classes == 'Lipid derivative') %>% pull(Gene_symbol)
vec_receptor_ncbi <- bitr(vec_receptor, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
ewiki <- enrichPathway(gene=vec_receptor_ncbi, pvalueCutoff = 0.05, readable=TRUE)
dim(ewiki)
head(ewiki, 3)

ewiki_head <- head(ewiki,10)
ewiki_head$Description <- factor(ewiki_head$Description, levels=rev(ewiki_head$Description))
ewiki_head$pvalue <- -log10(ewiki_head$p.adjust)

ggplot(ewiki_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Lipid derived HRs (Reactome Pathway)',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure1/pathway-reactome-Lipid-derivative.pdf', width = 6, height = 3)


# Peptide
vec_receptor <- df_receptor %>% filter(Hormone_chemical_classes == 'Peptide') %>% pull(Gene_symbol)
vec_receptor_ncbi <- bitr(vec_receptor, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
ewiki <- enrichPathway(gene=vec_receptor_ncbi, pvalueCutoff = 0.05, readable=TRUE)
dim(ewiki)
head(ewiki, 3)

ewiki_head <- head(ewiki,10)
ewiki_head$Description <- factor(ewiki_head$Description, levels=rev(ewiki_head$Description))
ewiki_head$pvalue <- -log10(ewiki_head$p.adjust)

ggplot(ewiki_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Peptide HRs (Reactome Pathway)',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure1/pathway-reactome-Peptide.pdf', width = 6, height = 3)


# Amino acid derivative
vec_receptor <- df_receptor %>% filter(Hormone_chemical_classes == 'Amino acid derivative') %>% pull(Gene_symbol)
vec_receptor_ncbi <- bitr(vec_receptor, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
ewiki <- enrichPathway(gene=vec_receptor_ncbi, pvalueCutoff = 0.05, readable=TRUE)
dim(ewiki)
head(ewiki, 3)

ewiki_head <- head(ewiki,10)
ewiki_head$Description <- factor(ewiki_head$Description, levels=rev(ewiki_head$Description))
ewiki_head$pvalue <- -log10(ewiki_head$p.adjust)

ggplot(ewiki_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Amino acid derived HRs (Reactome Pathway)',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure1/pathway-reactome-Amino-acid-derivative.pdf', width = 6, height = 3)


# orphan HR
vec_receptor <- df_receptor %>% filter(Hormone_chemical_classes == 'NA') %>% pull(Gene_symbol)
vec_receptor_ncbi <- bitr(vec_receptor, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
ewiki <- enrichPathway(gene=vec_receptor_ncbi, pvalueCutoff = 0.05, readable=TRUE)
dim(ewiki)
head(ewiki, 3)

ewiki_head <- head(ewiki,10)
ewiki_head$Description <- factor(ewiki_head$Description, levels=rev(ewiki_head$Description))
ewiki_head$pvalue <- -log10(ewiki_head$p.adjust)

ggplot(ewiki_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Orphan HRs (Reactome Pathway)',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
ggsave('output/HPA-single_cell/figure1/pathway-reactome-orphan-HR.pdf', width = 6, height = 3)

#################################
# Disease Ontology
vec_receptor_ncbi <- bitr(vec_receptor, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
go <- enrichDO(gene = vec_receptor_ncbi, ont = 'DO')
dim(go)

go_head <- head(go,20)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)

ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='Disease Ontology',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))

# ggsave('output/HPA-single_cell/figure1/pathway-DO.pdf', width = 4, height = 5)





#################################
# BP for various HR types classified by the corresponding hormone chemical classes
df_receptor$Hormone_chemical_classes <- recode(df_receptor$Hormone_chemical_classes, 'NA'='Orphan')
# Peptide
vec_receptor <- df_receptor %>% filter(Hormone_chemical_classes == 'Peptide') %>% pull(Gene_symbol)
go <- enrichGO(gene = vec_receptor, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
dim(go)
go_head <- head(go,5)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)
ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='GO biological process',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
# ggsave('output/HPA-single_cell/figure1/pathway-BP.pdf', width = 8, height = 5)

# Amino acid derivative
vec_receptor <- df_receptor %>% filter(Hormone_chemical_classes == 'Amino acid derivative') %>% pull(Gene_symbol)
go <- enrichGO(gene = vec_receptor, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
dim(go)
go_head <- head(go,5)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)
ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='GO biological process',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
# ggsave('output/HPA-single_cell/figure1/pathway-BP.pdf', width = 8, height = 5)

# Lipid derivative
vec_receptor <- df_receptor %>% filter(Hormone_chemical_classes == 'Lipid derivative') %>% pull(Gene_symbol)
go <- enrichGO(gene = vec_receptor, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
dim(go)
go_head <- head(go,5)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)
ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='GO biological process',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
# ggsave('output/HPA-single_cell/figure1/pathway-BP.pdf', width = 8, height = 5)

# Orphan
vec_receptor <- df_receptor %>% filter(Hormone_chemical_classes == 'Orphan') %>% pull(Gene_symbol)
go <- enrichGO(gene = vec_receptor, OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = 'BP')
dim(go)
go_head <- head(go,5)
go_head$Description <- factor(go_head$Description, levels=rev(go_head$Description))
go_head$pvalue <- -log10(go_head$p.adjust)
ggplot(go_head, aes(x=Count, y = Description, fill = pvalue)) +
  geom_bar(stat='identity', width = 0.8, alpha = 0.9) +
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 60),position = 'top')+
  geom_text(aes(label = Description), hjust = 0, size = 4, alpha = 1, x = 0.4)+
  #coord_flip() +
  scale_fill_gradient(low = '#E3E3E3',high = '#F9C6B0') +
  labs(x='Number of genes', y=NULL, fill = expression(paste(-log[10],italic(' p value')))) +
  guides(fill = guide_legend(title.position = "left"))+
  labs(title='GO biological process',hjust=0.5) +
  theme_classic()+
  scale_x_continuous(expand = c(0.01,0))+ #使Y轴与图形之间的缝隙消失
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.title=element_text(hjust=0.5,face="bold"),
        legend.position="left",
        legend.title = element_text(angle = 90, vjust= 0.9))
# ggsave('output/HPA-single_cell/figure1/pathway-BP.pdf', width = 8, height = 5)