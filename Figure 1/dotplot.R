library(Seurat)
library(ggplot2)
library(dplyr)
library(magrittr)
library(ggsci)

rm(list = ls())

data <- readRDS('seurat.RDS') # Only progenitors and bm cells
data@meta.data[which(data$condition == 'MNS'), ]$condition <- 'MNG' # MNG = multinodal goiter (English), MNS = multinodulair struma (Dutch)

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
