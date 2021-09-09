library(ggplot2)
library(openxlsx)
library(ggpubr)
library(dplyr)
library(magrittr)
library(rstatix)

setwd('/Users/martijnzoodsma/Documents/PhD/thyroid_cancer/AUC/')
rm(list = ls())
dev.off()

pbmc.bulk <- read.csv('enrichment_bulk_PBMC.csv')
bmmc.bulk <- read.csv('enrichment_bulk_BMMC.csv')

theme_set(theme_classic() + theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), plot.title = element_text(hjust = 0.5)))
celltypes.show <- c('B', 'CD14+ Monocyte', 'Naïve CD4+ T', 'CD8+ T', 'NK')
width = 0.8


# 1. Enrichment of BMMC DE genes in PBMC
names <- getSheetNames('enr_PBMC_cross.xlsx')

df2 <- data.frame()
for (name in names) { df2 <- rbind(df2, read.xlsx('enr_PBMC_cross.xlsx', sheet = name))}
colnames(df2) <- c('Tissue', 'BroadCelltype', 'Condition', 'AUC')
df2 %<>% filter(BroadCelltype %in% celltypes.show)
df2$BroadCelltype <- factor(df2$BroadCelltype, levels = celltypes.show)

# Filter out the MNS patients
df2 %<>% filter(Condition != 'MNS')

# TC greater than any of the others?
sig2 <- df2 %>% 
  group_by(Tissue, BroadCelltype) %>% 
  wilcox_test(AUC ~ Condition, alternative = 'greater', ref.group = 'TC') %>% 
  filter(p < 0.05) %>% 
  add_x_position(x = "BroadCelltype", group = 'Condition') %>% 
  add_y_position(step.increase = 0.1) %>% 
  add_significance(p.col = "p", output.col = 'label')

# Reorder x axis
df2$BroadCelltype <- factor(df2$BroadCelltype, levels = celltypes.show)

pl2 <- ggplot(data = df2) +
  geom_violin(aes(x = BroadCelltype, y = AUC, fill= Condition), width = width) +
  # geom_violin(data = bmmc.bulk, aes(x = 'Bulk', y = x.genes, fill = condition)) +
  ylim(c(0, 1)) +
  labs(title= 'Enrichment of BMMC DE genes in PBMC') +
  scale_x_discrete(drop=F) +
  theme(axis.title.x = element_blank()) +
  stat_pvalue_manual(sig2, label = 'label', tip.length = 0.01) +
  scale_fill_manual(values=c("#00AFBB", 
                             "#E7B800", 
                             "#FC4E07"))

pdf('../output/AUC/PBMC_enr.pdf', width = 6, height = 5, onefile=F)
pl2
dev.off()



# 2. Enrichment of PBMC DE genes in BMMC
names <- getSheetNames('enr_BMMC_cross.xlsx')

df2 <- data.frame()
for (name in names) { df2 <- rbind(df2, read.xlsx('enr_BMMC_cross.xlsx', sheet = name))}
colnames(df2) <- c('Tissue', 'BroadCelltype', 'Condition', 'AUC')

df2 %<>% filter(BroadCelltype %in% celltypes.show)
df2$BroadCelltype <- factor(df2$BroadCelltype, levels = celltypes.show)

df2 %<>% filter(Condition != 'MNS')


# Are TC higher enriched for cross-tissue TC-related differentially expressed genes?
sig3 <- df2 %>% 
  group_by(Tissue, BroadCelltype) %>% 
  wilcox_test(AUC ~ Condition, alternative = 'greater', ref.group = 'TC') %>% 
  filter(p < 0.05) %>% 
  add_x_position(x = "BroadCelltype", group = 'Condition') %>% 
  add_y_position(step.increase = 0.1) %>% 
  add_significance(p.col = "p", output.col = 'label')

pl3 <- ggplot(data = df2) +
  geom_violin(aes(x = BroadCelltype, y = AUC, fill = Condition), width = width) +
  # geom_violin(data = pbmc.bulk, aes(x = 'Bulk', y = x.genes, fill = condition), width = width) +
  ylim(c(0, 1)) +
  labs(title = 'Enrichment of PBMC DE genes in BMMC') +
  scale_x_discrete(drop=F) +
  theme(axis.title.x = element_blank()) +
  stat_pvalue_manual(sig3, label = 'label', tip.length = 0.01) +
  scale_fill_manual(values=c("#00AFBB", 
                             "#E7B800", 
                             "#FC4E07"))

pdf('../output/AUC/BMMC_enr.pdf', width = 6, height = 5, onefile=F)
pl3
dev.off()