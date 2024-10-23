#################################
# Author: Xiangjie Zhao
# Date: 2023.03.06
# Title: Analysis of hormone receptor gene list
#################################
rm(list = ls())
setwd('~/project/hormone_receptor_scRNA-2022-12-11/')

# Load Data
library("ggplot2") # for the plot
library("ggrepel") # for spreading text labels on the plot, you can replace with `geom_text` if you want
library("scales") # for axis labels notation
library('ggpubr')

# insert your steps to load data from tabular files or other sources here; 
# dummy datasets taken directly from files shown in this example
df_receptor <- readxl::read_excel('input/hormone_receptor_list.xlsx')

# hg38 chromosome sizes
chrom_sizes <- structure(list(Chromosome = rev(c("chr1", "chr2", "chr3", "chr4", 
                                             "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                                             "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                             "chr20", "chr21", "chr22", "chrX", "chrY")),
                              size = rev(c(249698942L, 242508799L, 198450956L, 190424264L, 181630948L,
                                       170805979L, 159345973L, 145138636L, 138688728L, 133797422L,
                                       135086622L, 133275309L, 114364328L, 108136338L, 102439437L,
                                       92211104L, 83836422L, 80373285L, 58617616L, 64444167L,
                                       46709983L, 51692466L, 156040895L, 57264655L))),
                         .Names = c("Chromosome", "size"),
                         class = "data.frame",
                         row.names = c(NA, -24L))

# hg38 centromere locations
# centromeres <- structure(list(Chromosome = c("chr1", "chr2", "chr3", "chr4", 
#                                              "chr5", "chr6", "chr7", "chr8", "chr9", "chrX", "chrY", "chr10", 
#                                              "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", 
#                                              "chr18", "chr19", "chr20", "chr21", "chr22"),
#                               start = c(121700000, 91800000, 87800000, 48200000, 46100000, 58500000,
#                                         58100000, 43200000, 42200000, 38000000, 51000000, 39254935L,
#                                         33200000, 16500000, 16100000, 17500000, 35300000, 22700000,
#                                         15400000, 24200000, 25700000, 10900000, 11288129L, 13000000L),
#                               end = c(125100000, 96000000, 94000000, 51800000, 51400000, 62600000,
#                                       62100000, 47200000, 45500000, 41600000, 55800000, 42254935L,
#                                       37800000, 18900000, 18200000, 20500000, 38400000, 27400000,
#                                       21500000, 28100000, 30400000, 29369569L, 14288129L, 16000000L)),
#                          .Names = c("Chromosome", "start", "end"),
#                          class = "data.frame",
#                          row.names = c(NA, -24L))
centromeres <- structure(list(Chromosome = rev(c('chr1','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr2','chr20','chr21','chr22','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrX','chrY')),
                              start = rev(c(121700000,38000000,51000000,33200000,16500000,16100000,17500000,35300000,22700000,15400000,24200000,91800000,25700000,10900000,13700000,87800000,48200000,46100000,58500000,58100000,43200000,42200000,58100000,10300000)),
                              end = rev(c(125100000,41600000,55800000,37800000,18900000,18200000,20500000,38400000,27400000,21500000,28100000,96000000,30400000,13000000,17400000,94000000,51800000,51400000,62600000,62100000,47200000,45500000,63800000,10600000))
                              ),
                         .Names = c("Chromosome", "start", "end"),
                         class = "data.frame",
                         row.names = c(NA, -24L))

#################################
# Adjust Data
# create an ordered factor level to use for the chromosomes in all the datasets
chrom_order1 <- rev(c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                 "chr22", "chrX", "chrY"))
chrom_key <- setNames(object = as.character(rev(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                                              12, 13, 14, 15, 16, 17, 18, 19, 20, 
                                              21, 22, 23, 24))), nm = chrom_order1)
chrom_order <- factor(x = chrom_order1, levels = chrom_order1)

# convert the chromosome column in each dataset to the ordered factor
chrom_sizes[["Chromosome"]] <- factor(x = chrom_sizes[["Chromosome"]], levels = chrom_order1)
df_receptor[["Chromosome"]] <- factor(x = df_receptor[["Chromosome"]], levels = chrom_order1)
centromeres[["Chromosome"]] <- factor(x = centromeres[["Chromosome"]], levels = chrom_order1)
# create a color key for the plot
group.colors <- c(Nucleus = "#237EB6", Membrane = "#ED1B28")

#################################
# Make plot
ggplot(data = chrom_sizes) + 
  # base rectangles for the chroms, with numeric value for each chrom on the x-axis
  geom_rect(aes(xmin = as.numeric(Chromosome) - 0.2, 
                xmax = as.numeric(Chromosome) + 0.2, 
                ymax = size, ymin = 0), 
            colour="black", fill = "white") + 
  # rotate the plot 90 degrees
  coord_flip() +
  # black & white color theme 
  theme(axis.text.x = element_text(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank()) + 
  # give the appearance of a discrete axis with chrom labels
  scale_x_discrete(name = "Chromosome", limits = names(chrom_key)) +
  # add bands for centromeres
  geom_rect(data = centromeres, aes(xmin = as.numeric(Chromosome) - 0.2, 
                                    xmax = as.numeric(Chromosome) + 0.2, 
                                    ymax = end, ymin = start)) +
  # add bands for Receptor_subcellular value
  geom_rect(data = df_receptor, aes(xmin = as.numeric(Chromosome) - 0.2, 
                                   xmax = as.numeric(Chromosome) + 0.2, 
                                   ymax = End+400000, ymin = Start-400000,
                                   fill = Receptor_subcellular), alpha=0.8) + 
  scale_fill_manual(values = group.colors) +
  # add 'Nucleus' gene markers
  geom_text_repel(data = subset(df_receptor, df_receptor$Receptor_subcellular == "Nucleus"), 
                  aes(x = Chromosome, y = Start, label = Gene_symbol), 
                  color = "#237EB6", show.legend = FALSE, size=1.7, max.overlaps = 200, alpha=0.9) +
  # add 'Membrane' gene markers
  geom_text_repel(data = subset(df_receptor, df_receptor$Receptor_subcellular == "Membrane"), 
                  aes(x = Chromosome, y = Start, label = Gene_symbol ), 
                  color = "#ED1B28", show.legend = FALSE, size=1.7, max.overlaps = 200, alpha=0.9) +
  ggtitle("Distribution of HRs on chromosomes") +
  # supress scientific notation on the y-axis
  scale_y_continuous(labels = comma) +
  ylab("Region (bp)") +
  labs(fill='HR types')

# ggsave('output/HPA-single_cell/figure1/chromosome-distribution1.pdf', width = 7, height = 6)



#################################
#统计染色体上的HR个数
library(tidyverse)
df <- df_receptor %>% group_by(Chromosome, Receptor_subcellular) %>% summarise(count = n()) %>% ungroup()
df %>% add_row(Chromosome='chrY', Receptor_subcellular='Membrane', count=0) %>% add_row(Chromosome='chrY', Receptor_subcellular='Nucleus', count=0) -> df
df %>% add_row(Chromosome='chr21', Receptor_subcellular='Membrane', count=0) %>% add_row(Chromosome='chr21', Receptor_subcellular='Nucleus', count=0) -> df
df$Chromosome <- factor(df$Chromosome, levels = chrom_order)
ggplot(df, aes(Chromosome, count, fill=Receptor_subcellular))+
  #geom_bar(position = 'stack') +
  geom_bar(stat = 'identity') +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c('#ED1B28','#237EB6')) +
  labs(y='Count')

# ggsave('output/HPA-single_cell/figure1/chromosome-distribution-count.pdf', width = 3, height = 5.75)

#统计染色体上的核/膜HR的密度
chrom_sizes_ocuppy <- df_receptor %>% group_by(Chromosome, Receptor_subcellular) %>% summarise(size_occupy = max(End)-min(Start)) %>% ungroup()

df_dens <- left_join(df, chrom_sizes_ocuppy, by=c('Chromosome', 'Receptor_subcellular'))
df_dens <- df_dens %>% filter(count > 1) # 对于某种类型的HR，某染色体上只有一个或0个基因时不计算密度，舍去
df_dens <- df_dens %>% mutate(dens = count/size_occupy)

mycoms <- list(c('Membrane', 'Nucleus'))
ggplot(df_dens, aes(x=Receptor_subcellular, y=log10(dens), fill=Receptor_subcellular)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.3, size=0.5)+
  stat_compare_means(method = 'wilcox.test', comparisons = mycoms)+
  theme_classic() +
  scale_fill_manual(values = c('#ED1B28','#237EB6')) +
  labs(x='', y='Density (log10)', fill='HR types')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
# ggsave('output/HPA-single_cell/figure1/chromosome-distribution-density.pdf', width = 2.5, height = 2.5)


#统计染色体上的核/膜HR的密度（用滑窗方法计算密度）
#xxxxxxx



#################################
#比较各个染色体上的激素受体密度
df <- df %>% group_by(Chromosome) %>% summarise(s = sum(count)) %>% ungroup()
df <- df[order(df$Chromosome),]
chrom_sizes <- chrom_sizes[order(chrom_sizes$Chromosome),]
chr_density_1e8 <- df$s / chrom_sizes$size * 1e8
df$chr_density_1e8 <- chr_density_1e8
# df$Chromosome <- factor(df$Chromosome, levels = chrom_order)
df <- df %>% arrange(desc(chr_density_1e8))
df$Chromosome <- factor(df$Chromosome, levels = df$Chromosome)

ggplot(df, aes(Chromosome, chr_density_1e8)) +
  geom_bar(stat = 'identity', width = 0.8, fill='white', color='black') +
  # coord_flip() +
  theme_classic() +
  labs(y='HR density in chromosome (No./100Mbps)') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"))

ggsave('output/HPA-single_cell/figure1/HR-density-in-chromosome-1.pdf', width = 9, height = 3.5)

