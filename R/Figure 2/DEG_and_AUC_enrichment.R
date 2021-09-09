library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(magrittr)
library(AUCell)
library(DESeq2)

# Cross-tissue comparison
  # 1. Calculate DE genes per celltype in PBMC and BMMC separately
  # 2. Calculate bulk DE genes (background)
  # 3. In PBMC:
    # Check enrichment of the bulk DE genes in pseudobulk from BMMC
    # Check enrichment of celltype specific DE genes per celltype from BMMC
  # 4. Idem dito in BMMC





rm(list = ls())
setwd('/vol/projects/mzoodsma/thyroid_cancer/analysis')
OUTDIR <- '/vol/projects/mzoodsma/thyroid_cancer/output/DEG/'


data <- readRDS('../data.RDS')
Idents(data) <- 'condition'

data@meta.data[which(data@meta.data$condition == 'MNS'), ]$condition <- 'MNG'
print(table(data$Tissue, data$condition))




# 1. Calculate DE genes per tissue for each celltype separately
# PBMC first

pbmc.counts <- matrix(NA, nrow = length(unique(data$BroadCelltype)), ncol = 2)
colnames(pbmc.counts) <- c('Downregulated', 'Upregulated')
rownames(pbmc.counts) <- unique(data$BroadCelltype)

for (ct in unique(data$BroadCelltype)) {
  
  # Calculate DE genes
  sub.pbmc <- subset(x = data, subset = BroadCelltype == ct & Tissue == 'PBMC')
  print(table(sub.pbmc$Tissue, sub.pbmc$BroadCelltype))
  
  de <- FindMarkers(object = sub.pbmc, ident.1 = 'TC', ident.2 = 'HC', test.use = 'MAST',
                    latent.vars = c('age', 'sex'), logfc.threshold = 0.0, min.cells.feature = 1, min.cells.group = 1)
  de$gene <- rownames(de)
  de.genes <- de %>% filter(p_val_adj < 0.05) %>% pull(gene)
  
  pbmc.counts[ct, 'Downregulated'] <- de %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC < 0) %>% nrow()
  pbmc.counts[ct, 'Upregulated'] <- de %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0) %>% nrow()
  
  write.csv(file = paste0(OUTDIR, 'PBMC/', ct, '_DE.csv'), de)
}
write.csv(file = paste0(OUTDIR, 'PBMC/DEG_counts.csv'), pbmc.counts)

# BMMC
bmmc.counts <- matrix(NA, nrow = length(unique(data$BroadCelltype)), ncol = 2)
colnames(bmmc.counts) <- c('Downregulated', 'Upregulated')
rownames(bmmc.counts) <- unique(data$BroadCelltype)

for (ct in unique(data$BroadCelltype)) {
  
  # Calculate DE genes
  sub.pbmc <- subset(x = data, subset = BroadCelltype == ct & Tissue == 'BMMC')
  print(table(sub.pbmc$Tissue, sub.pbmc$BroadCelltype))
  
  de <- FindMarkers(object = sub.pbmc, ident.1 = 'TC', ident.2 = 'HC', test.use = 'MAST',
                    latent.vars = c('age', 'sex'), logfc.threshold = 0.0, min.cells.feature = 1, min.cells.group = 1)
  de$gene <- rownames(de)
  de.genes <- de %>% filter(p_val_adj < 0.05) %>% pull(gene)
  
  bmmc.counts[ct, 'Downregulated'] <- de %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC < 0) %>% nrow()
  bmmc.counts[ct, 'Upregulated'] <- de %>% filter(p_val_adj < 0.05) %>% filter(avg_log2FC > 0) %>% nrow()
  
  write.csv(file = paste0(OUTDIR, 'BMMC/', ct, '_DE.csv'), de)
}
write.csv(file = paste0(OUTDIR, 'BMMC/DEG_counts.csv'), bmmc.counts)


# Calculate DE genes for bulk for each tissue separately
df <- matrix(NA, nrow = nrow(data@assays$RNA@counts),
             ncol = length(unique(data$sample)),
             dimnames = list(rownames(data@assays$RNA@counts),
                             unique(data$sample))
)

for (s in unique(data$sample)) {
  cells <- data@meta.data %>% filter(sample == s) %>% rownames(.)
  counts <- data@assays$RNA@counts[ , cells]
  expr <- Matrix::rowSums(counts)
  df[, s] <- expr
}

meta <- data@meta.data %>% select(all_of(c('age', 'condition', 'sample', 'sex', 'Tissue'))) %>% distinct()
print(table(meta$condition, meta$Tissue))
meta$sex <- factor(meta$sex)
meta$Tissue <- factor(meta$Tissue)

# Calculate DE genes for PBMC and BMMC separately
# Include PCA plots for Tissue, sex and conditon
pbmc.meta <- meta %>% filter(Tissue == 'PBMC')
rownames(pbmc.meta) <- pbmc.meta$sample
pbmc.df <- df[, pbmc.meta$sample]
stopifnot(all(colnames(pbmc.df) == rownames(pbmc.meta)))

dds<- DESeqDataSetFromMatrix(countData = pbmc.df,
                             colData = pbmc.meta,
                             design = ~ condition + sex)


vsd <- varianceStabilizingTransformation(dds)
theme_set(theme_classic() + theme(plot.title = element_text(hjust = 0.5)))

pdf(paste0(OUTDIR, 'PBMC/PCA_bulk.pdf'), width = 10, height = 3)
annotate_figure(
  ggarrange(
    plotPCA(vsd, intgroup = 'condition') + labs(title ='Condition'),
    plotPCA(vsd, intgroup = 'sex') + labs(title ='sex'),
    plotPCA(vsd, intgroup = 'Tissue') + labs(title ='Tissue'), nrow = 1),
  top = 'PBMC pseudobulk')
dev.off()

pdf(paste0(OUTDIR, 'PBMC/PCA_bulk_condition.pdf'), width = 5, height = 5)
plotPCA(vsd, intgroup = 'condition') + labs(title = 'PBMC')
dev.off()
stop()

dds <- DESeq(dds)
TC_vs_HC <- results(dds, contrast = c('condition', 'TC', 'HC')) %>% as.data.frame() %>% arrange(padj)
TC_vs_MNG <- results(dds, contrast = c('condition', 'TC', 'MNG')) %>% as.data.frame() %>% arrange(padj)
MNG_vs_HC <- results(dds, contrast = c('condition', 'MNG', 'HC')) %>% as.data.frame() %>% arrange(padj)

write.csv(file = paste0(OUTDIR, 'PBMC/bulk_DE_TC_vs_HC.csv'), TC_vs_HC)
write.csv(file = paste0(OUTDIR, 'PBMC/bulk_DE_TC_vs_MNG.csv'), TC_vs_MNG)
write.csv(file = paste0(OUTDIR, 'PBMC/bulk_DE_MNG_vs_HC.csv'), MNG_vs_HC)


# BMMC pseudobulk genes
bmmc.meta <- meta %>% filter(Tissue == 'BMMC')
rownames(bmmc.meta) <- bmmc.meta$sample
bmmc.df <- df[, bmmc.meta$sample]
stopifnot(all(colnames(bmmc.df) == rownames(bmmc.meta)))

dds<- DESeqDataSetFromMatrix(countData = bmmc.df,
                             colData = bmmc.meta,
                             design = ~ condition + sex)
vsd <- varianceStabilizingTransformation(dds)
theme_set(theme_classic() + theme(plot.title = element_text(hjust = 0.5)))

pdf(paste0(OUTDIR, 'BMMC/PCA_bulk.pdf'), width = 10, height = 3)
annotate_figure(
  ggarrange(
    plotPCA(vsd, intgroup = 'condition') + labs(title ='Condition'),
    plotPCA(vsd, intgroup = 'sex') + labs(title ='sex'),
    plotPCA(vsd, intgroup = 'Tissue') + labs(title ='Tissue'), nrow = 1),
  top = 'BMMC pseudobulk')
dev.off()

pdf(paste0(OUTDIR, 'BMMC/PCA_bulk_condition.pdf'), width = 4, height = 4)
plotPCA(vsd, intgroup = 'condition') + labs(title = 'BMMC')
dev.off()

dds <- DESeq(dds)
TC_vs_HC <- results(dds, contrast = c('condition', 'TC', 'HC')) %>% as.data.frame() %>% arrange(padj)
TC_vs_MNG <- results(dds, contrast = c('condition', 'TC', 'MNG')) %>% as.data.frame() %>% arrange(padj)
MNG_vs_HC <- results(dds, contrast = c('condition', 'MNG', 'HC')) %>% as.data.frame() %>% arrange(padj)

write.csv(file = paste0(OUTDIR, 'BMMC/bulk_DE_TC_vs_HC.csv'), TC_vs_HC)
write.csv(file = paste0(OUTDIR, 'BMMC/bulk_DE_TC_vs_MNG.csv'), TC_vs_MNG)
write.csv(file = paste0(OUTDIR, 'BMMC/bulk_DE_MNG_vs_HC.csv'), MNG_vs_HC)

stop()
OUTDIR <- '/vol/projects/mzoodsma/thyroid_cancer/output/AUC/'
# AUC enrichment
# 1. Check the enrichment of the cross-tissue pseudobulk genes
# 2. Cell-type specific: Check genes from cross-tissue and same tissue as check

# PBMC
# Enrichment of BMMC pseudobulk genes
df <- read.csv(paste0('/vol/projects/mzoodsma/thyroid_cancer/output/DEG/BMMC/bulk_DE.csv'), row.names=1)
de.genes.bulk <- df %>% filter(padj < 0.05) %>% rownames(.)

# Enrichment of bulk DE from BMMC in bulk PBMC
sub.pbmc <- subset(x = data, subset = Tissue == 'PBMC')
expr <- GetAssayData(sub.pbmc)
geneset = list(genes = de.genes.bulk)
rankings <- AUCell_buildRankings(exprMat = expr, nCores = 5, plotStats = F)
cells.auc <- AUCell_calcAUC(geneSets = geneset, rankings = rankings, nCores = 5, aucMaxRank = nrow(rankings) * 0.05)
auc <- cells.auc@assays@data@listData$AUC
x <- t(cells.auc@assays@data@listData$AUC) %>% data.frame()
stopifnot(all(rownames(x) == rownames(sub.pbmc@meta.data)))
exp <- cbind(sub.pbmc@meta.data %>% select(all_of(c('Tissue', 'BroadCelltype', 'condition'))), x$genes)
write.csv(file = paste0(OUTDIR, 'PBMC/enrichment_bulk_BMMC.csv'), exp)

results.pbmc.self <- list()
results.pbmc.cross <- list()

for (ct in unique(data$BroadCelltype)) {
  pbmc.genes <- read.csv(paste0('/vol/projects/mzoodsma/thyroid_cancer/output/DEG/PBMC/', ct, '_DE.csv')) %>% filter(p_val_adj < 0.05) %>% pull(gene)
  bmmc.genes <- read.csv(paste0('/vol/projects/mzoodsma/thyroid_cancer/output/DEG/BMMC/', ct, '_DE.csv')) %>% filter(p_val_adj < 0.05) %>% pull(gene)
  sub.pbmc <- subset(x = data, subset = BroadCelltype == ct & Tissue == 'PBMC')
  expr <- GetAssayData(sub.pbmc)
  
  # Enrichment of celltype-specific DE from PBMC: Self
  if (length(pbmc.genes > 0)) {
    geneset = list(genes = pbmc.genes)
    rankings <- AUCell_buildRankings(exprMat = expr, nCores = 5, plotStats = F)
    cells.auc <- AUCell_calcAUC(geneSets = geneset, rankings = rankings, nCores = 5, aucMaxRank = nrow(rankings) * 0.05)
    auc <- cells.auc@assays@data@listData$AUC
    x <- t(cells.auc@assays@data@listData$AUC) %>% data.frame()
    stopifnot(all(rownames(x) == rownames(sub.pbmc@meta.data)))
    exp <- cbind(sub.pbmc@meta.data %>% select(all_of(c('Tissue', 'BroadCelltype', 'condition'))), x$genes)
    print(head(exp))
    results.pbmc.self[[ct]] <- exp
  }
  
  # Enrichment of celltype-specific DE from BMMC: Cross
  if (length(bmmc.genes > 0)) {
    geneset = list(genes = bmmc.genes)
    rankings <- AUCell_buildRankings(exprMat = expr, nCores = 5, plotStats = F)
    cells.auc <- AUCell_calcAUC(geneSets = geneset, rankings = rankings, nCores = 5, aucMaxRank = nrow(rankings) * 0.05)
    auc <- cells.auc@assays@data@listData$AUC
    x <- t(cells.auc@assays@data@listData$AUC) %>% data.frame()
    stopifnot(all(rownames(x) == rownames(sub.pbmc@meta.data)))
    
    exp <- cbind(sub.pbmc@meta.data %>% select(all_of(c('Tissue', 'BroadCelltype', 'condition'))), x$genes)
    print(head(exp))
    results.pbmc.cross[[ct]] <- exp
  }
}
openxlsx::write.xlsx(x = results.pbmc.cross, file = paste0(OUTDIR, 'PBMC/enr_PBMC_cross.xlsx'))
openxlsx::write.xlsx(x = results.pbmc.self, file = paste0(OUTDIR, 'PBMC/enr_PBMC_self.xlsx'))


# BMMC
# Enrichment of PBMC pseudobulk genes in PBMC
df <- read.csv(paste0('/vol/projects/mzoodsma/thyroid_cancer/output/DEG/PBMC/bulk_DE.csv'), row.names=1)
de.genes.bulk <- df %>% filter(padj < 0.05) %>% rownames(.)

sub.bmmc <- subset(x = data, subset = Tissue == 'BMMC')
expr <- GetAssayData(sub.bmmc)
geneset = list(genes = de.genes.bulk)
rankings <- AUCell_buildRankings(exprMat = expr, nCores = 5, plotStats = F)
cells.auc <- AUCell_calcAUC(geneSets = geneset, rankings = rankings, nCores = 5, aucMaxRank = nrow(rankings) * 0.05)
auc <- cells.auc@assays@data@listData$AUC
x <- t(cells.auc@assays@data@listData$AUC) %>% data.frame()
stopifnot(all(rownames(x) == rownames(sub.bmmc@meta.data)))
exp <- cbind(sub.bmmc@meta.data %>% select(all_of(c('Tissue', 'BroadCelltype', 'condition'))), x$genes)
write.csv(file = paste0(OUTDIR, 'BMMC/enrichment_bulk_PBMC.csv'), exp)

results.bmmc.self <- list()
results.bmmc.cross <- list()

for (ct in unique(data$BroadCelltype)) {
  pbmc.genes <- read.csv(paste0('/vol/projects/mzoodsma/thyroid_cancer/output/DEG/PBMC/', ct, '_DE.csv')) %>% filter(p_val_adj < 0.05) %>% pull(gene)
  bmmc.genes <- read.csv(paste0('/vol/projects/mzoodsma/thyroid_cancer/output/DEG/BMMC/', ct, '_DE.csv')) %>% filter(p_val_adj < 0.05) %>% pull(gene)
  sub.bmmc <- subset(x = data, subset = BroadCelltype == ct & Tissue == 'BMMC')
  expr <- GetAssayData(sub.bmmc)
  
  if (length(pbmc.genes > 0)) {
    # Enrichment of celltype-specific DE from PBMC: Cross
    geneset = list(genes = pbmc.genes)
    rankings <- AUCell_buildRankings(exprMat = expr, nCores = 5, plotStats = F)
    cells.auc <- AUCell_calcAUC(geneSets = geneset, rankings = rankings, nCores = 5, aucMaxRank = nrow(rankings) * 0.05)
    auc <- cells.auc@assays@data@listData$AUC
    x <- t(cells.auc@assays@data@listData$AUC) %>% data.frame()
    stopifnot(all(rownames(x) == rownames(sub.bmmc@meta.data)))
    exp <- cbind(sub.bmmc@meta.data %>% select(all_of(c('Tissue', 'BroadCelltype', 'condition'))), x$genes)
    
    results.bmmc.cross[[ct]] <- exp
  }
  
  if (length(bmmc.genes > 0)) {
    # Enrichment of celltype-specific DE from BMMC: Self
    geneset = list(genes = bmmc.genes)
    rankings <- AUCell_buildRankings(exprMat = expr, nCores = 5, plotStats = F)
    cells.auc <- AUCell_calcAUC(geneSets = geneset, rankings = rankings, nCores = 5, aucMaxRank = nrow(rankings) * 0.05)
    auc <- cells.auc@assays@data@listData$AUC
    x <- t(cells.auc@assays@data@listData$AUC) %>% data.frame()
    stopifnot(all(rownames(x) == rownames(sub.bmmc@meta.data)))
    exp <- cbind(sub.bmmc@meta.data %>% select(all_of(c('Tissue', 'BroadCelltype', 'condition'))), x$genes)
    results.bmmc.self[[ct]] <- exp
    
  }
}

openxlsx::write.xlsx(x = results.bmmc.cross, file = paste0(OUTDIR, 'BMMC/enr_BMMC_cross.xlsx'))
openxlsx::write.xlsx(x = results.bmmc.self, file = paste0(OUTDIR, 'BMMC/enr_BMMC_self.xlsx'))
