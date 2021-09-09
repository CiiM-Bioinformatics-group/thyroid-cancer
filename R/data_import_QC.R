# Loading libraries
library(Seurat)
library(ggplot2)
library(plyr)

# Project title"
OUTPUT_DIR=paste0("/Users/kuijpers/Documents/PhD/SepsisThyroid/TC/plots/Both/pipeline/")

sampID <- list()
sampID[[4]] <- c("0" = "TC18.02", "1" = "TC18.08", "2" = "TC18.21", "3" = "TC18.52", "4" = "TC18.42")
sampID[[5]] <- c("0" = "TC18.36", "1" = "TC18.41", "2" = "TC18.01", "3" = "TC18.51", "4" = "TC18.26", "5" = "TC18.07")
sampID[[6]] <- c("0" = "TC18.25", "1" = "TC18.37", "2" = "TC18.55", "3" = "TC18.45", "4" = "TC18.23", "5" = "TC18.01")
sampID[[7]] <- c("0" = "TC18.12", "1" = "TC18.11", "2" = "TC18.38", "3" = "TC18.22", "4" = "TC18.23", "5" = "TC18.04")
sampID[[8]] <- c("0" = "TC18.31", "1" = "TC18.07", "2" = "TC18.24", "3" = "TC18.11", "4" = "TC18.54")
sampID[[9]] <- c("0" = "TC18.46", "1" = "TC18.22", "2" = "TC18.02", "3" = "TC18.26", "4" = "TC18.35")
sampID[[10]] <- c("0" = "TC18.06", "1" = "TC18.24", "2" = "TC18.43", "3" = "TC18.03", "4" = "TC18.57", "5" = "TC18.25")

obj <- list()
qcData <- data.frame(Pool=character(6), all=character(6), dplx=character(6), qc=character(6), stringsAsFactors = F)
pools <- c(4:10)
for (q in 1:length(pools)) {
  p <- pools[q]
  # Reading in the cluster assignment
  clusters <- read.table(paste0("/Users/kuijpers/Documents/PhD/SepsisThyroid/TC/Data/Pool",p,"/", "clusters.tsv"), header = TRUE, stringsAsFactors = FALSE)
  clusters$barcode <- substr(clusters$barcode, start = 1, stop = nchar(clusters$barcode)-2)
  rownames(clusters) <- clusters$barcode
  doublets <- subset(clusters, status == "doublet")
  unassigned <- subset(clusters, status == "unassigned")
  clusters <- subset(clusters, status == "singlet")
  # Reading in the 10X data
  seur <- Read10X(paste0("/Users/kuijpers/Documents/PhD/SepsisThyroid/TC/Data/Pool",p,"/filtered_feature_bc_matrix/"))
  # We remove all cells that are not singlets
  colnames(seur) <- substring(colnames(seur),1,nchar(colnames(seur))-2)
  seur <- seur[, colnames(seur) %in% clusters$barcode]
  # We calculate the percentage of MT reads
  MT.index <- grep(pattern = "^MT-", x = rownames(seur), value = FALSE)
  all.sum <- Matrix::colSums(seur)
  percent.MT <- Matrix::colSums(seur[MT.index, ])/all.sum
  seur <- seur[-MT.index, ] # remove MT genes
  # We calculate the percentage of Ribosomal genes
  RG.index <- grep("^RP[SL][0-9]+$", x = rownames(seur), perl=TRUE, value = FALSE)
  percent.RG <- Matrix::colSums(seur[RG.index, ])/all.sum
  seur <- seur[-RG.index, ] # remove Ribosomal genes
  # Adding the cluster assignment
  clusters <- clusters[,c(1:3)]
  clusters$sample <- as.factor(clusters$assignment)
  sampMatch <- sampID[[p]]
  clusters$sample <- as.character(revalue(clusters$sample, sampMatch))
  cluster.assignment <- clusters$sample
  names(cluster.assignment) <- clusters$barcode
  cluster.assignment <- cluster.assignment[order(match(names(cluster.assignment), names(percent.MT)))]
  # Adding the status assignment	
  cluster.status <- clusters$status
  names(cluster.status) <- clusters$barcode
  cluster.status <- cluster.status[order(match(names(cluster.status), names(percent.MT)))]
  # Adding the meta data
  seur <- CreateSeuratObject(seur, min.cells = 5, meta.data = data.frame(percent.mt = percent.MT, percent.rg = percent.RG, sample = cluster.assignment, status = cluster.status), project = paste0("Pool",p))
  seur <- subset(seur, nFeature_RNA > 250 & nFeature_RNA < 6000 & percent.mt < 0.2)
  seur$pool <- p
  seur$Tissue <- "PBMC"
  
  obj[[paste0("POOL",p)]] <- seur
}

# merge items in obj then delete obj
data <- merge(obj[["POOL4"]], y = c(obj[["POOL5"]], obj[["POOL6"]], obj[["POOL7"]], obj[["POOL8"]], obj[["POOL9"]], obj[["POOL10"]]), add.cell.ids = c("Pool4", "Pool5", "Pool6","Pool7", "Pool8", "Pool9", "Pool10"), project = "Thyroid")
rm(sampID, sampMatch, obj, clusters, doublets, qcData, seur, unassigned, all.sum, allCells, cluster.assignment, cluster.status, dplexCells, MT.index, p, percent.MT, percent.RG, pools, q, qcCells, RG.index)

data@meta.data[which(data$sample %in% c("TC18.52") & data$pool == 4),]$Tissue <- "BMMC"
data@meta.data[which(data$sample %in% c("TC18.51", "TC18.01") & data$pool == 5),]$Tissue <- "BMMC"
data@meta.data[which(data$sample %in% c("TC18.23", "TC18.55") & data$pool == 6),]$Tissue <- "BMMC"
data@meta.data[which(data$sample %in% c("TC18.11", "TC18.12") & data$pool == 7),]$Tissue <- "BMMC"
data@meta.data[which(data$sample %in% c("TC18.07","TC18.24","TC18.54") & data$pool == 8),]$Tissue <- "BMMC"
data@meta.data[which(data$sample %in% c("TC18.02","TC18.26","TC18.22") & data$pool == 9),]$Tissue <- "BMMC"
data@meta.data[which(data$sample %in% c("TC18.25","TC18.57") & data$pool == 10),]$Tissue <- "BMMC"

# Visualize QC metrics as a violin plot
Idents(data) <- data$pool
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rg"), ncol = 2)

###########################################################################################
#add some basic data
YLINKED_GENES="FAM197Y5,TSPY6P,TSPY10,AC016752.1,CDY1B,BPY2B,FAM197Y1,SLC9B1P1,DAZ3,USP9Y,DAZ4,DDX3Y,UTY,BPY2C,CDY1,AC007965.1,AC010877.1,TMSB4Y,VCY,VCY1B,NLGN4Y,CDY2B,CDY2A,HSFY1,HSFY2,AC009977.1,KDM5D,EIF1AY,RPS4Y2,CYorf17,RBMY1B,RBMY1A1,RBMY1E,RBMY1D,SRY,RPS4Y1,PRY2,RBMY1F,ZFY,TGIF2LY,PCDH11Y,AC012067.1,TSPY2,RBMY1J,PRY,BPY2,AMELY,DAZ1,TBL1Y,TSPY4,TSPY8,TSPY3,TSPY1,DAZ2"
ylinked <- unlist(strsplit(YLINKED_GENES, ","))
XIC="XIST"
xic <- unlist(strsplit(XIC, ","))

sampInfo <- read.table("/Users/kuijpers/Documents/PhD/SepsisThyroid/TC/Data/TC.patient.info.txt", sep = "\t", stringsAsFactors = F, header = T)
rownames(sampInfo) <- sampInfo$ID

data@meta.data$age <- NA
data@meta.data$sex <- NA
data@meta.data$condition <- NA
data@meta.data$diagnosis <- NA
data@meta.data$treatment <- NA

for (s in sampInfo$ID) {
  try({
    data@meta.data[which(data@meta.data$sample == s),]$age <- sampInfo[s,2]
    data@meta.data[which(data@meta.data$sample == s),]$sex <- sampInfo[s,3]
    data@meta.data[which(data@meta.data$sample == s),]$condition <- sampInfo[s,4]
    data@meta.data[which(data@meta.data$sample == s),]$diagnosis <- sampInfo[s,5]
    data@meta.data[which(data@meta.data$sample == s),]$treatment <- sampInfo[s,6]
  })
}

Male <- subset(data, sex %in% c("m"))
Male$sex <- "Male"
Female <- subset(data, sex %in% c("f"))
Female$sex <- "Female"

# We determine which female cells have expressed y-linked genes
Female.ylinked <- colSums(as.data.frame(Female[ylinked, ]@assays$RNA@counts, stringsAsFactors = FALSE))
# Now, we remove the y-linked genes epressing cells from our female samples
Female.ylinked <- Female.ylinked[Female.ylinked == 0]
# Now we subset the female cells to only include the cells not expressing ylinked genes
Female <- subset(Female, cells = names(Female.ylinked))

# We determine which male cels have expressed the Xic genes
Male.xic <- colSums(as.data.frame(Male[xic, ]@assays$RNA@counts, stringsAsFactors = FALSE))
# Now, we remove the y-linked genes epressing cells from our female samples
Male.xic <- Male.xic[Male.xic == 0]
# Now we subset the female cells to only include the cells not expressing ylinked genes
Male <- subset(Male, cells = names(Male.xic))

# Now, we merge the male and female cells
data <- merge(x = Male, y = Female)
rm(Female, Male)

###########################################################################################


# Plotting the expression of y-linked genes
pdf(paste0(OUTPUT_DIR, "dist/HeatMap_Y-linked.pdf"), width = 12, height = 7)
print(DoHeatmap(data, features = ylinked, group.by = "sample", slot = "data") +
        theme(axis.text = element_text(size = 20)) +
        ggtitle(label = "Expression of y-linked genes") +
        scale_fill_gradient2(low = "#4575B4", mid = "#FEFEC0", high = "#D73027", midpoint = 0))
dev.off()


###########################################################################################

# Processing data for the UMAP
#data <- SCTransform(data, vars.to.regress = c("percent.mt", "percent.rg"))

data <- NormalizeData(data) #norm
data <- ScaleData(data, vars.to.regress = c("percent.mt", "percent.rg"))
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- RunPCA(data, npcs = 20, verbose = FALSE)

# Generate the UMAP of the raw Seurat Object
data <- RunUMAP(data, reduction = "pca", dims = 1:20)
data <- FindNeighbors(data, reduction = "pca", dims = 1:20)
data <- FindClusters(data, resolution = 0.65, algorithm = 3, n.start = 20, n.iter = 20)

###########################################################################################


# Plotting the UMAP of the different samples
pdf(paste0(OUTPUT_DIR, "dist/UMAP_samples.pdf"), width = 12, height = 7)
print(DimPlot(data, group.by = "sample", pt.size = 0.5, reduction = "umap") +
        ggtitle(label = "UMAP by Samples") +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()

pdf(paste0(OUTPUT_DIR, "dist/UMAP_pools.pdf"), width = 12, height = 7)
print(DimPlot(data, group.by = "pool", pt.size = 0.5, reduction = "umap") +
        ggtitle(label = "UMAP by batch") +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()

# UMAP cannot illustrate the difference between the conditions
# We therefore plot a PCA
pdf(paste0(OUTPUT_DIR, "dist/PCA_samples.pdf"), width = 12, height = 7)
print(DimPlot(data, group.by = "sample", pt.size = 3, reduction = "pca") +
        ggtitle(label = "PCA by Samples") +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()

pdf(paste0(OUTPUT_DIR, "dist/PCA_pools.pdf"), width = 12, height = 7)
print(DimPlot(data, group.by = "pool", pt.size = 3, reduction = "pca") +
        ggtitle(label = "PCA by batch") +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()

###########################################################################################

# Plotting the Seurat clustering
pdf(paste0(OUTPUT_DIR, "dist/UMAP_SeuratClusters.pdf"), width = 12, height = 7)
print(DimPlot(data, group.by = "seurat_clusters", pt.size = 0.5, reduction = "umap", label = T) + NoLegend() +
        ggtitle(label = "UMAP of Seurat Clusters") +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()

pdf(paste0(OUTPUT_DIR, "dist/PCA_SeuratClusters.pdf"), width = 12, height = 7)
print(DimPlot(data, group.by = "seurat_clusters", pt.size = 3, reduction = "pca") +
        ggtitle(label = "PCA by Samples") +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()

###########################################################################################

# Plotting the cells by condition
pdf(paste0(OUTPUT_DIR, "dist/UMAP_Condition.pdf"), width = 12, height = 7)
print(DimPlot(data, group.by = "condition", pt.size = 0.5, reduction = "umap",
              cols = c("green","lightgrey", "red")) +
        ggtitle(label = "UMAP of Condition") +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()

pdf(paste0(OUTPUT_DIR, "dist/PCA_Condition.pdf"), width = 12, height = 7)
print(DimPlot(data, group.by = "condition", pt.size = 3, reduction = "pca",
              cols = c("green","lightgrey", "red")) +
        ggtitle(label = "PCA by Condition") +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()

###########################################################################################

# We save the results of the QC so we can use it in other scripts
save.image(paste0("/Users/kuijpers/Documents/PhD/SepsisThyroid/TC/plots/Both/data/Thyroid.RData"))

