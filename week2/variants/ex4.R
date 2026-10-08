#Load libraries
library(ggplot2)
library(tidyverse)

#Load data
gt_long <- read.delim("/Users/cmdb/qb26-answers/week2/variants/gt_long.txt", header=FALSE, sep= "\t")
AF <- read.delim("/Users/cmdb/qb26-answers/week2/variants/AF.txt", header=FALSE, sep ="\t")

#Label the columns more nicely
colnames(AF) <- c("Chr", "Pos", "Frequency")

#plot
ggplot(AF) +
  geom_histogram(aes(x = Frequency), bins = 11) +
  xlab("Allele Frequency") +
  ylab("Count")

setwd("/Users/cmdb/qb26-answers/week2/variants")
ggsave("AF.png")

colnames(gt_long) <- c("chromosome", "position", "sample", "genotype")

#ChrII of sample A01_62, 
A01_62_ChrII <- gt_long %>% filter(chromosome == "chrII", sample == "A01_62")
ggplot(A01_62_ChrII, aes(x = position, y = sample, color = genotype)) +
  geom_point() +
  xlab("Position") +
  ylab("Genotype")

#All chromosomes on A01_62
A01_62_all <- gt_long %>%
  filter(sample == "A01_62")
ggplot(A01_62_all, aes(x = position, y = sample, color = genotype)) +
  geom_point() +
  facet_grid("chromosome")
xlab("Position") +
  ylab("Genotype")
  
  
#All chromosomes on all samples
ggplot(gt_long, aes(x = position, y = sample, color = genotype)) +
  geom_point() +
  facet_grid(.~chromosome) +
  xlab("Position") +
  ylab("Genotype")

ggsave("ancestry.png", width=45, height=15 )
