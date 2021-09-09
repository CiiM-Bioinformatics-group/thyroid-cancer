## GO enrichment in CD14+ Monocytes, Naive CD4+ T cells and CD8+ T cells

# Fig. 2C
rm(list = ls())

library(ggplot2)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(openxlsx)

# Heatmap of the DEg in CD14+, CD8+ T and Naive CD4
setwd('/Users/martijnzoodsma/Documents/PhD/thyroid_cancer/')

# CD14+ Monocytes
df1 <- read.csv('rsync_output_server/DEG/PBMC/CD14+ Monocyte/TC_vs_HC_DEG.csv') %>% filter(p_val_adj < 0.05) %>% filter(abs(avg_log2FC)>0.5) %>% mutate(celltype = 'CD14+ Monocyte') %>% mutate(comparison = 'TC vs HC')
df2 <- read.csv('rsync_output_server/DEG/PBMC/CD14+ Monocyte/MNG_vs_HC_DEG.csv') %>%  filter(p_val_adj < 0.05) %>% filter(abs(avg_log2FC)>0.5)%>% mutate(celltype = 'CD14+ Monocyte') %>% mutate(comparison = 'MNG vs HC')
df3 <- read.csv('rsync_output_server/DEG/PBMC/CD14+ Monocyte/TC_vs_MNG_DEG.csv') %>%  filter(p_val_adj < 0.05) %>% filter(abs(avg_log2FC)>0.5)%>% mutate(celltype = 'CD14+ Monocyte') %>% mutate(comparison = 'TC vs MNG')

# CD8+ T cells
df4 <- read.csv('rsync_output_server/DEG/PBMC/CD8+ T/TC_vs_HC_DEG.csv') %>% filter(p_val_adj < 0.05) %>% filter(abs(avg_log2FC)>0.5) %>% mutate(celltype = 'CD8+ T') %>% mutate(comparison = 'TC vs HC')
df5 <- read.csv('rsync_output_server/DEG/PBMC/CD8+ T/MNG_vs_HC_DEG.csv') %>% filter(p_val_adj < 0.05) %>% filter(abs(avg_log2FC)>0.5) %>% mutate(celltype = 'CD8+ T') %>% mutate(comparison = 'MNG vs HC')
df6 <- read.csv('rsync_output_server/DEG/PBMC/CD8+ T/TC_vs_MNG_DEG.csv') %>% filter(p_val_adj < 0.05) %>% filter(abs(avg_log2FC)>0.5) %>% mutate(celltype = 'CD8+ T') %>% mutate(comparison = 'TC vs MNG')

# Naive CD4+ T cells
df7 <- read.csv('rsync_output_server/DEG/PBMC/Naïve CD4+ T/TC_vs_HC_DEG.csv') %>% filter(p_val_adj < 0.05) %>% filter(abs(avg_log2FC)>0.5) %>% mutate(celltype = 'Naïve CD4+ T') %>% mutate(comparison = 'TC vs HC')
df8 <- read.csv('rsync_output_server/DEG/PBMC/Naïve CD4+ T/MNG_vs_HC_DEG.csv') %>% filter(p_val_adj < 0.05) %>% filter(abs(avg_log2FC)>0.5) %>% mutate(celltype = 'Naïve CD4+ T') %>% mutate(comparison = 'MNG vs HC')
df9 <- read.csv('rsync_output_server/DEG/PBMC/Naïve CD4+ T/TC_vs_MNG_DEG.csv') %>% filter(p_val_adj < 0.05) %>% filter(abs(avg_log2FC)>0.5) %>% mutate(celltype = 'Naïve CD4+ T') %>% mutate(comparison = 'TC vs MNG')

# All genes
####
# This is not used any more
####

# genes <- unique(c(df1$gene, df2$gene, df3$gene, df4$gene, df5$gene, df6$gene, df7$gene, df8$gene, df9$gene))
# 
# df <- rbind(df1 ,df2, df3, df4, df5, df6, df7, df8, df9) %>% select(gene, avg_log2FC, celltype, comparison) 
# df$gene <- factor(df$gene, levels = rev(unique(df$gene)))
# df$comparison = factor(df$comparison, levels = c('TC vs HC', 'MNG vs HC', 'TC vs MNG'))
# df <- reshape2::melt(df)
# df$xlabel <- paste0(df$celltype, ' ', df$comparison)
# 
# pdf('output/fig 2C.pdf', height = 8, width = 5)
# ggplot(data = df) +
#   geom_tile(aes(x = comparison, y = gene, fill = value), width = 0.5) +
#   theme_classic() +
#   theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) +
#   facet_wrap(~celltype) +
#   labs(x = 'Comparison', y = 'Gene', fill = 'log2FC') +scale_fill_gradient(low = 'darkblue', high = 'red')
# dev.off()


#############################
# Enrichment of the celltypes Fig 2D plot only

rm(list = ls())
dev.off()
setwd('/Users/martijnzoodsma/Documents/PhD/thyroid_cancer/rsync_output_server/DEG/PBMC/')

res <- list(
  'CD14+ Mono_MNG'      = read.xlsx('CD14+ Monocyte/MNG_vs_HC_enr.xlsx', sheet = 'GO upreg') %>% filter(Count >= 5), 
  'CD14+ Mono_TC'       = read.xlsx('CD14+ Monocyte/TC_vs_HC_enr.xlsx', sheet = 'GO upreg')%>% filter(Count >= 5),
  'CD8+ T_MNG'          = read.xlsx('CD8+ T/MNG_vs_HC_enr.xlsx', sheet = 'GO upreg') %>% filter(Count >= 5),
  'CD8+ T_TC'           = read.xlsx('CD8+ T/TC_vs_HC_enr.xlsx', sheet = 'GO upreg')%>% filter(Count >= 5),
  'Naive CD4+ T_MNG'    = read.xlsx('Naïve CD4+ T/MNG_vs_HC_enr.xlsx', sheet = 'GO upreg') %>% filter(Count >= 5),
  'Naive CD4+ T_TC'     = read.xlsx('Naïve CD4+ T/TC_vs_HC_enr.xlsx', sheet = 'GO upreg')%>% filter(Count >= 5)
) %>% merge_result()

df <- res@compareClusterResult
df$old_cluster <- df$Cluster
ids <- df$old_cluster
df %<>% select(-Cluster)

df <- cbind(
  df, 
  reshape2::colsplit(string = ids, pattern = '_', names = c('Cluster', 'Condition'))
)

df$GeneRatio2 <- sapply(df$GeneRatio, function(x) eval(parse(text=x)))
df <- df %>% group_by(old_cluster) %>% slice(1:5)
df$Description <- factor(df$Description, levels = unique(df$Description))
df %<>% arrange(Cluster, Condition)

pdf('../../../output/enrichment_celltypes_separate_TC_MNG.pdf', width = 10, height = 6)
ggplot(data = df) +
  geom_point(aes(x = Cluster, y = Description, size = GeneRatio2, color = p.adjust)) +
  DOSE::theme_dose() +
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), 
        strip.background = element_rect(fill = 'white', color = 'black'), axis.title.y = element_blank(), axis.title.x = element_blank()) +
  scale_colour_gradient(limits=c(0, 0.05), low="red", high = 'blue') +
  facet_wrap('Condition')
dev.off()

openxlsx::write.xlsx('GOenr_celltypes_conditions_TC.xlsx', x = df)

# No longer used. We use GO enrichment only to reduce the number of subfigures
# 
# # KEGG enrichment plot
# df1 <- read.xlsx('CD14+ Monocyte/MNG_vs_HC_enr.xlsx', sheet = 'KEGG up') %>% filter(Count >= 5)
# df2 <- read.xlsx('CD14+ Monocyte/TC_vs_HC_enr.xlsx', sheet = 'KEGG up')%>% filter(Count >= 5)
# df3 <- read.xlsx('CD8+ T/TC_vs_HC_enr.xlsx', sheet = 'KEGG up')%>% filter(Count >= 5)
# df4 <- read.xlsx('Naïve CD4+ T/MNG_vs_HC_enr.xlsx', sheet = 'KEGG up')%>% filter(Count >= 5)
# 
# res <- list(
#   'CD14+ Mono MNG' = df1, 
#   'CD14+ Mono TC' = df2, 
#   'CD8+ T' = df3, 
#   'Naive CD4+ T' = df4
# )
# 
# pdf('../../../output/fig_2D_KEGG.pdf', width = 7, height = 6)
# dotplot(merge_result(res)) +
#   theme(axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45), plot.title = element_text(hjust = .5)) +
#   labs(title = 'Upregulated KEGG pathways')
# dev.off()
