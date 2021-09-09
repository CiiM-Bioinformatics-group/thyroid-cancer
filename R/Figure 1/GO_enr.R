rm(list = ls())

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(magrittr)
library(openxlsx)
library(Seurat)

enrich_all <- function(up, down, outpath_xlsx, outpath_pdf) {
  
  enr.GO.up <- enrichGO(gene = up,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 0.05,
                        readable = TRUE)
  
  enr.GO.down <- enrichGO(gene = down,
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 1,
                          qvalueCutoff = 0.05,
                          readable = TRUE)
  
  enr.KEGG.up <- enrichKEGG(gene = up,
                            organism = 'hsa',
                            pAdjustMethod = "BH",
                            pvalueCutoff = 1,
                            qvalueCutoff = 0.05)
  
  enr.KEGG.down <- enrichKEGG(gene = down,
                              organism = 'hsa',
                              pAdjustMethod = "BH",
                              pvalueCutoff = 1,
                              qvalueCutoff = 0.05)
  
  res = list('GO upreg' = enr.GO.up,
             'GO downreg' = enr.GO.down,
             'KEGG up' = enr.KEGG.up,
             'KEGG down' = enr.KEGG.down)
  write.xlsx(x = res, file = outpath_xlsx)
  
  pdf(outpath_pdf, width = 15, height = 15)
  print(dotplot(merge_result(res), showCategory = 10) +
          theme(text = element_text(size =12),
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                plot.title = element_text(hjust = 0.5)))
  dev.off()
  
}


data <- readRDS('../data.RDS')
data <- subset(data, Tissue == 'PBMC')
data@meta.data[which(data$condition == 'MNS'), ]$condition <- 'MNG'

nr_DEGs <- matrix(0, nrow = length(unique(data@meta.data$BroadCelltype)), ncol = 2)
colnames(nr_DEGs) = c('Downregulated', 'Upregulated')
rownames(nr_DEGs) = unique(data@meta.data$BroadCelltype)
Idents(data) <- 'condition'
OUTPUT_DIR <- '/vol/projects/mzoodsma/thyroid_cancer/output/DEG/celltypes/'

total <- list()

for (cell in unique(data@meta.data$BroadCelltype)) {
  
  sub <- subset(data, BroadCelltype == cell)
  
  # Use wilcoxon test for now. MAST when age and gender of the samples are known
  # Set ident.1 as after vaccination, and ident.2 as before vaccination. This way,
  # logFC values > 0 indicate upregulated genes AFTER vaccination.
  res <- FindMarkers(object = sub, ident.1 = 'TC', ident.2 = 'HC', test.use = 'MAST', latent.vars = c('age', 'sex'), logfc.threshold = 0.0, min.cells.feature = 1, min.cells.group = 1)
  
  res$significance <- ifelse(res$p_val_adj < 0.05, TRUE, FALSE)
  res$direction <- ifelse(res$avg_log2FC > 0, 'Upregulated', 'Downregulated')
  res %>% arrange(p_val_adj)
  res$gene <- rownames(res)
  
  # Export
  write.csv(res, paste0(OUTPUT_DIR, cell, '_DEG.csv'))
  
  # Count number of sig genes and export this as well
  sig.res <- res %>% filter(p_val_adj < 0.05)
  nr_DEGs[cell, 'Downregulated'] <- sig.res %>% filter(avg_log2FC < 0) %>% nrow()
  nr_DEGs[cell, 'Upregulated'] <- sig.res %>% filter(avg_log2FC > 0) %>% nrow()
  
  # Volcano plot per celltype
  # Top X significant labeled
  pdf(paste0(OUTPUT_DIR, cell, '_volcano.pdf'), width = 5, height = 5)
  print(ggplot() +
          geom_point(data = res %>% filter(significance == F) , aes(x = avg_log2FC, y = -log10(p_val_adj)), color = 'lightgray', show.legend = F) +
          geom_point(data = res %>% filter(significance == T) , aes(x = avg_log2FC, y = -log10(p_val_adj), color = direction)) +
          labs(x = 'Log-fold change after / before vaccination', y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = paste0(cell)) +
          scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
          geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
          geom_vline(xintercept = 0.0, lty = 'dotted') +
          theme_classic() +
          theme(plot.title = element_text(hjust = 0.5)) +
          ggrepel::geom_label_repel(data = head(x = res, 15), aes(label = gene, x = avg_log2FC, y = -log10(p_val_adj)), box.padding = 0.5, show.legend = F),
        
  )
  dev.off()
  
  # Enrichment
  enrich_all(up = sig.res %>% filter(p_val_adj < 0.05) %>% filter(direction == 'Upregulated') %>% pull(gene),
             down = sig.res %>% filter(p_val_adj < 0.05) %>% filter(direction == 'Downregulated') %>% pull(gene),
             outpath_xlsx = paste0(OUTPUT_DIR, 'enr_', cell, '.xlsx'),
             outpath_pdf = paste0(OUTPUT_DIR, 'volcano_', cell, '.pdf'))
}

write.csv(paste0(OUTPUT_DIR, 'DEG/nr_DEGs.csv'), x = nr_DEGs)
