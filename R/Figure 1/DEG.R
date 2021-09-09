library(Seurat)
library(ggplot2)
library(dplyr)
library(magrittr)
library(ggsci)
library(EnhancedVolcano)
library(ggpubr)

rm(list = ls())

data <- readRDS('seurat.RDS') # Only progenitors and bm cells
data@meta.data[which(data$condition == 'MNS'), ]$condition <- 'MNG' # MNG = multinodal goiter (English), MNS = multinodulair struma (Dutch)

# DE analyses
# Within multipotent progenitors. All conditions vs all
unique(data2$celltype_new)
data3 <- subset(data2, celltype_new == 'Myeloid/multipotent progenitors')
Idents(data3) <- 'condition'

x1 <- FindMarkers(object = data3, ident.1 = 'TC', ident.2 = 'HC', test.use = 'MAST', latent.vars = c('age', 'sex'), logfc.threshold = 0.0, min.cells.feature = 1, min.cells.group = 1)
x2 <- FindMarkers(object = data3, ident.1 = 'TC', ident.2 = 'MNG', test.use = 'MAST', latent.vars = c('age', 'sex'), logfc.threshold = 0.0, min.cells.feature = 1, min.cells.group = 1)
x3 <- FindMarkers(object = data3, ident.1 = 'MNG', ident.2 = 'HC', test.use = 'MAST', latent.vars = c('age', 'sex'), logfc.threshold = 0.0, min.cells.feature = 1, min.cells.group = 1)

# Nr of DEGS per comparison
x1 %>% filter(p_val_adj < 0.05) %>% nrow()
x2 %>% filter(p_val_adj < 0.05) %>% nrow()
x3 %>% filter(p_val_adj < 0.05) %>% nrow()

pdf('../output/DE_MPP.pdf', width = 10 ,height = 5, onefile = F)
ggarrange(
  EnhancedVolcano(toptable = x1, x = 'avg_log2FC', y = 'p_val_adj', lab = rownames(x1), 
                  pCutoff = 0.05, subtitle = 'TC vs HC', title = element_blank()),
  EnhancedVolcano(toptable = x2, x = 'avg_log2FC', y = 'p_val_adj', lab = rownames(x2), 
                  pCutoff = 0.05, subtitle = 'MNG vs HC', title = element_blank()),
  EnhancedVolcano(toptable = x3, x = 'avg_log2FC', y = 'p_val_adj', lab = rownames(x3), 
                  pCutoff = 0.05, subtitle = 'TC vs MNG', title = element_blank()), 
  nrow = 1, common.legend = T
)

dev.off()
