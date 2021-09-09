library(Seurat)
library(ggplot2)
library(dplyr)
library(magrittr)
library(AUCell)
library(ggsci)

rm(list = ls())

data <- readRDS('seurat.RDS') # Only progenitors and bm cells
data@meta.data[which(data$condition == 'MNS'), ]$condition <- 'MNG' # MNG = multinodal goiter (English), MNS = multinodulair struma (Dutch)

###############################################################################
# UMAPs coloured by several factors
pdf('../output/UMAP_condition.pdf', width = 12, height = 7)
print(DimPlot(data, group.by = 'condition', pt.size = 0.5, reduction = 'umap') +
        ggtitle('umap by condition') +
        theme(text = element_text(size = 15),
              axis.text = element_text(size = 10)))
dev.off()

pdf('../output/UMAP_tissue.pdf', width = 12, height = 7)
umap_condition <- DimPlot(data, group.by = 'Tissue', pt.size = 0.5, reduction = 'umap') +
  ggtitle('umap by condition') +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 10))

print(umap_condition)
dev.off()

pdf('../output/UMAP_celltype.pdf', width = 12, height = 7)
print(DimPlot(data, group.by = 'celltype', pt.size = 0.5, reduction = 'umap') +
        ggtitle('umap by celltype') +
        theme(text = element_text(size = 15),
              axis.text = element_text(size = 10)))
dev.off()

pdf('../output/UMAP_broadcelltype.pdf', width = 12, height = 7)
print(DimPlot(data, group.by = 'BroadCelltype', pt.size = 0.5, reduction = 'umap') +
        ggtitle('umap by celltype') +
        theme(text = element_text(size = 15),
              axis.text = element_text(size = 10)))
dev.off()

pdf('../output/UMAP_sub_celltype.pdf', width = 12, height = 7)
print(DimPlot(data, group.by = 'sub.celltypes', pt.size = 0.5, reduction = 'umap') +
        ggtitle('umap by celltype') +
        theme(text = element_text(size = 15),
              axis.text = element_text(size = 10)))
dev.off()

#### Rename the celltype column to the same names used in the plots
data@meta.data$celltype_new <- data@meta.data$celltype
data@meta.data[which(data@meta.data$celltype == 'CLP'), 'celltype_new'] <- 'Lymphoid progenitors'
data@meta.data[which(data@meta.data$celltype == 'MPP'), 'celltype_new'] <- 'Myeloid/multipotent progenitors'
data@meta.data[which(data@meta.data$celltype == 'unclassified'), 'celltype_new'] <- 'Early B cell progenitors'
data@meta.data[which(data@meta.data$celltype == 'Pre-B Cell'), 'celltype_new'] <- 'Late B cell progenitors'

pdf('../output/UMAP_celltype_new.pdf', width = 12, height = 7)
umap_celltype <- DimPlot(data, group.by = 'celltype_new', pt.size = 0.5, reduction = 'umap') +
  ggtitle('umap by celltype') +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 10))

print(umap_celltype)
dev.off()



###############################################################################
# Dotplot of the GWAS-identified TC-related genes
data2 <- subset(data, celltype_new != 'B - PBMC')
print(table(data2$celltype_new))

data2$celltype_new <- factor(data2$celltype_new, levels = rev(c('Late B cell progenitors', 'Early B cell progenitors', 'Lymphoid progenitors', 'Myeloid/multipotent progenitors')))
Idents(data2) <- 'celltype_new'

print(table(data2$celltype))
print(table(data2$celltype_new))

features <- c('FLT3', 'KIT', 'SPN', 'CD34', 'MME', 'IL7R', 'CD24', 'CD38', 'CD19', 'PTPRC', 'CD79A', 'MS4A1', 'IGHM', 'IGHD', 'NCAM1', 'CD14', 'NCR1', 'ITGAX', 'IL4R', 'CR2', 'CD40')

pdf('../output/dotplot_literature_genes.pdf', width = 10, height = 5)
print(DotPlot(data2,features=features, cols = "RdBu") +
        theme(legend.position = "top",legend.title = element_blank()) +
        xlab("") +
        ylab("") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        ggtitle(label = "Bone marrow derived celltype markers"))
dev.off()




