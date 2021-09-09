# This script annotates the demultiplexed cells
# This step requires a lot of user input 



# Reading in the filtered Seurat object
load(paste0("/Users/kuijpers/Documents/PhD/SepsisThyroid/TC/plots/Both/data/Thyroid.RData"))


# Loading the necessary libraries
library(Seurat)
library(ggplot2)
library(gridExtra)

##################################################################################################

# Violin Plot of MT and ERCC intensity
get.violin.data1 <- function(seurat, genes) {
  #  output = data.frame(gene = character(0), value= numeric(0), ident = character(0))
  #####################################################################################################
  output = data.frame(gene = character(0), value= numeric(0), ident = character(0), tech=character(0))
  #####################################################################################################
  for (gene in genes) {
    if(any(gene == seurat@assays$RNA@data@Dimnames[[1]])){
      data.use = data.frame(FetchData(seurat,gene))
      data.use = t(data.use)
      data.melt=data.frame(rep(gene, length(seurat@active.ident)))
      colnames(data.melt)[1]="gene"
      data.melt$value=as.numeric(data.use[1,1:length(seurat@active.ident)])
      data.melt$id=names(data.use)[1:length(seurat@active.ident)]
      data.melt$ident=seurat@active.ident
      ############################################
      #data.melt$tech=seurat$stim # ???What was used to stimulate the cells?
      ############################################
      if(any(data.melt$value != 0)) noise = rnorm(length(data.melt$value))/100000 else noise = 0
      data.melt$value=as.numeric(as.character(data.melt$value))+noise
      output = rbind(output, data.melt)
    } else {
      data.melt=data.frame(
        gene = rep(gene, seurat@assays$RNA@data@Dim[2]),
        value = rep(0, seurat@assays$RNA@data@Dim[2]), 
        ident = seurat@active.ident
      )
      output = rbind(output, data.melt)
    }
  }
  return(output)
}

##################################################################################################

# Marker Genes
# We use the following list of marker genes to annotate our cells
# Source: Suppressive myeloid cells are a hallmark of severe COVID-19, Cell 2020
# IL7R,CCR7,TCF7  Naive CD4 T
# IL7R            Memory CD4 T
# CD8A,GZMK       CD8 T
# NKG7,GZMB       NK
# CD14, LYZ       CD14+ monocytes
# FCGR3A          CD16+ monocytes
# CST3, CD86      MDC, DC
# IL1B            proinflammatory
# CD79A           B cell
# CD27,SDC1       plasma
# PPBP            Mk
# KIT,TPSAB1      Mast

#PBMC narkers
# Source: obmc3k vigenette seurat
# Markers 	      Cell Type
# IL7R, CCR7 	    Naive CD4+ T
# IL7R, S100A4 	  Memory CD4+
# CD14, LYZ 	    CD14+ Mono
# MS4A1 	        B
# CD8A 	          CD8+ T
# FCGR3A, MS4A7 	FCGR3A+ Mono
# GNLY, NKG7 	    NK
# FCER1A, CST3 	  DC
# PPBP 	          Platelet

#some gene markers for general PBMC annotation,
#IL7R+,CCR7+,TCF7+                            Naive CD4 T  
#IL7R,CCR7-,TCF7-                             Memory CD4 T   
#CD8A+,CD8B+,GZMK+                            CD8 T    
#NKG7+,NCAM1+,CD8A-                           NK    
#CD14+, LYZ+, FCGR3A-                         CD14+ monocytes   
#FCGR3A+, CD14-                               CD16+ monocytes   
#CST3+, CD86+, HLA-DRA+, CD14-, FCGR3A-       mDC, pDC (dendritic cells)    
#IL1B,TNF,IL6                                 proinflammatory markers  # macrophages
#CD79A+, CD27-, SDC1-                         B cell        
#CD79A+, CD27+, SDC1+                         plasmablast        
#PPBP+                                        megakaryocyte/platelet   
#KIT+, TPSAB1+                                Mast

# Pre-Neutrophiles, Immature Neutrophiles and Mature Neutrophiles all express the following genes
# S100A8, S100A9
# Pre-Neutrophiles express the following genes
# ELANE, MPO, PRTN3
# Pre-Neutrophiles and Immature Neutrophiles express the following genes
# PADI4, CD24, ANXA1, BPI, LCN2, ARG1, GSN, OLFM4
# Mature Neutrophiles do not express most of the marker genes. They do express S100A8 and S100A9.
# Mature Neutrophiles express the following genes, but only at low levels.
# PADI4, ANXA1, ARG1, GSN
# CD274 is only lowly expressed by mature neutrophiles.
# Source: Figure 4, Suppressive myeloid cells are a hallmark of severe COVID-19, Cell 2020


#MNC Population	Markers expressed (positive) or not expressed (negative)								
#                                                  CD3 
#                                                  CD19
#                                                  CD56
#                                                  CD10
#                                                  CD14
#                                                  CD66b
#                                                  CD335

#                                                  CD11c	  CD34	      CD38	   CD45RA   	CD90    	CD10	    CD123   	  CD110   	CD135
#HSC  Hematopoietic Stem Cell	                    negative  positive	negative	negative	positive	negative      -           -         -
#MPP	Multipotent Progenitor	                    negative  positive	negative	negative	negative	negative      -           -         -
#LMPP	Lymphoid-Primed Multipotent Progenitor	    negative  positive	negative	positive	negative	negative      -           -         -
#MLP	Multi Lymphoid Progenitor	                  negative  positive	negative	positive	negative	positive	    -           -         -	

#CMP	Common Myeloid Progenitor	                  negative  positive	positive	negative	negative	negative	intermediate	  -	      positive
#MEP	Megakaryocyte-Erythrocyte Progenitor	      negative  positive	positive	negative	negative	negative	negative	    positive	negative
#GMP	Granulocyte-Macrophage Progenitor	          negative	positive  positive	positive	negative	negative	intermediate	  -	      positive
#CLP	Common Lymphoid Progenitor	                negative  positive	positive	positive	negative	positive	    -           -         -

# celltype_marker <- c("IL7R","CCR7","TCF7","CD8A","GZMK","NKG7","GZMB","CD14","LYZ","FCGR3A","CST3","CD86","IL1B","CD79A","CD27","SDC1","PPBP","KIT","TPSAB1",
#                      "S100A8","S100A9","ELANE","MPO","PRTN3","PADI4","CD24","ANXA1","BPI","LCN2","ARG1","GSN","OLFM4", "CD274")
markers_pbmc <- c("CD34","CD45","IL7R","CCR7","S100A4","TCF7","CD8A","GZMK","GNLY","NKG7","GZMB","CD14","LYZ","MS4A1","MS4A7","FCER1A","FCGR3A","CST3","CD86","IL1B",
                  "CD79A","CD27","SDC1","PPBP","KIT","TPSAB1")
markers_neutrophil <- c("S100A8","S100A9","ELANE","MPO","PRTN3","PADI4","CD24","ANXA1","BPI","LCN2","ARG1","GSN","OLFM4", "CD274")
markers_bmmc <- c("CD45", "CD19", "CD10", "CD20", "CD27", "CD21", "CD38", "CD138", "CD11c", "CD123", "CD3", "CD14", "CD16","HLA-DR", "CD141", "CD1c", "CD34", "CD117", 
                  "CD56", "CD33", "CD13", "CD11B", "CD64", "CD7", "CD8", "CD45RA", "CD69", "CD103", "CD4")
markers_mnc <- c("CD3", "CD19", "CD56", "CD10", "CD14", "CD66b", "CD335", "CD11c", "CD34", "CD38", "CD45RA", "CD90", "CD123", "CD110", "CD135")

markers_covid19 <- c("CCR7", "LEHF1", "FHIT", "PIK3IP1", "LDHB", "MAL", "ADTRP", "EEF1B2", "RPS3A", "TCF7", "LTB", "IL7R", "AQP3", "JUNB", "TRADD", "RPSA", "RPLP0", 
                     "RPS6", #CD4 T 18 markers
                     "VCAN", "LYZ", "S100A9", "S100A12", "LGALS1", "MAFB", "CD14", "PLBD1", "FCN1", "S100A8", #Classical Monocytes 10 markers (total until here 27)
                     "GZMH", "CCL5", "NKG7", "FGFBP2", "CD8A", "GZMA", "GNLY", "CST7", "CD8B", "CCL4", # CD8 T 10 markers (total until here 38)
                     "GZMB", "PRF1", "SPON2", "CLIC3", "CTSW", "HOPX", # NK cells 6 markers (total until here 44)
                     "IFII27", "IFITM3", "APOBEC3A", "TYMP", "ISG15", "CST3", "LGALS2", "IFI6", # Inflammatory Monocytes 8 markers (total until here 52)
                     "CD79A", "MS4A1", "CD74", "CD79B", "HLA-DRA", "TCL1A", "HLA-DQB1", "HLA-DPA1", "HLA-DQA1", "HLA-DPB1", # B cells 10 markers (total until here 62)
                     "GZMK", "DUSP2", "LYAR", "KLRG1", "IL32", "KLRB1", "S100B", "RPS5", "NELL2", "NOSIP", # CD 8 T 10 markers (total until here 72)
                     "FCGR3B", "NAMPT", "IFITM2", "CXCL8", "NEAT1", "G0S2", "SLC25A37", "PROK2", "BASP1", "CSF3R", # Neutrophils 10 markers (total until here 82)
                     "CORO1B", "CD27", "RTKN2", "ARID5B", "TRAC", "ITGB1", "TRBC2", "TRBC1", "IFI44L", "XAF1", "MX1", "IFIT3", "STAT1", "EIF2AK2", # CD4 T 14 markers (total until here 96)
                     "LST1", "FCGR3A", "AIF1", "MS4A7", "C1QA", "COTL1", "CFD", "PSAP", "FTL", # Non-classical Monocytes 9 markers (total until here 105)
                     "LTF", "LCN2", "CAMP", "RETN", "DEFA4", "CD24", "PGLYRP1", "MMP8", # Immature Neutrophils 8 markers (total until here 113)
                     "AC02O656.1", "TNFAIP2", "MPEG1", "CPVL", "FGL2", "CSTA", "CYBB", "MTRNR2L12", "CTSS", # Monocytes (ZFP36L2) 9 markers (total until here 122)
                     "STMN1", "MKI67", "TYMS", "PCNA", "TUBA1B", "TUBB", "DUT", "HMBG2", "HMGN2", "HIST1H4C", # Prol. T cells 10 markers (total until here 132)
                     "FCER1A", "HLA-DRB5", "HLA-BRB1", # mDCs 3 markers (total until here 135)
                     "JCHAIN", "IGHG1", "IGHG3", "MYB1", "IGHA1", "IGLC7", "IGLC2", "IGHM", "IGLC3", "IGKC", # Plasmablast 10 markers (total until here 145)
                     "PPBP", "PF4", "NRGN", "GNG11", "TUBB1", "CAVIN2", "HBG2", "MYL9", "GP9", "CLU", # Megakaryocytes 10 markers (total until here 155)
                     "HBA2", "HBA1", "HBD", "ALAS2", "HBM", "AHSP", "CA1", "SNCA", "HBB", # Erythroblasts 9 markers (total until here 164)
                     "SOX4", "PRSS57", "SPINK2", "AVP", "EGFL7", "CDK6", "ANKRD28", "SERPINB1", "AREG", # Undefined 9 markers (total until here 173)
                     "ITM2C", "PLD4", "SERPINF1", "TCF4", "LILRA4", "IRF8", "IRF7", "PTGDS", "PPP1R14B", # pDCs 9 markers (total until here 182)
                     "CMTM2", "SOD2", "PLEK" # CD8 3 markers (total until here 185)
)
marker.list <- list(
  c("CCR7", "LEHF1", "FHIT", "PIK3IP1", "LDHB", "MAL", "ADTRP", "EEF1B2", "RPS3A", "TCF7", "LTB", "IL7R", "AQP3", "JUNB", "TRADD", "RPSA", "RPLP0", "RPS6"),
  c("CORO1B", "CD27", "RTKN2", "ARID5B", "TRAC", "ITGB1", "TRBC2", "TRBC1", "IFI44L", "XAF1", "MX1", "IFIT3", "STAT1", "EIF2AK2"),
  c("VCAN", "LYZ", "S100A9", "S100A12", "LGALS1", "MAFB", "CD14", "PLBD1", "FCN1", "S100A8"),
  c("GZMH", "CCL5", "NKG7", "FGFBP2", "CD8A", "GZMA", "GNLY", "CST7", "CD8B", "CCL4", "GZMK", "DUSP2", "LYAR", "KLRG1", "IL32", "KLRB1", "CMTM2", "SOD2", "PLEK"),
  c("GZMB", "PRF1", "SPON2", "CLIC3", "CTSW", "HOPX"),
  c("IFII27", "IFITM3", "APOBEC3A", "TYMP", "ISG15", "CST3", "LGALS2", "IFI6"),
  c("CD79A", "MS4A1", "CD74", "CD79B", "HLA-DRA", "TCL1A", "HLA-DQB1", "HLA-DPA1", "HLA-DQA1", "HLA-DPB1"),
  c("S100B", "RPS5", "NELL2", "NOSIP"),
  c("FCGR3B", "NAMPT", "IFITM2", "CXCL8", "NEAT1", "G0S2", "SLC25A37", "PROK2", "BASP1", "CSF3R"),
  c("LST1", "FCGR3A", "AIF1", "MS4A7", "C1QA", "COTL1", "CFD", "PSAP", "FTL"),
  c("LTF", "LCN2", "CAMP", "RETN", "DEFA4", "CD24", "PGLYRP1", "MMP8"),
  c("AC02O656.1", "TNFAIP2", "MPEG1", "CPVL", "FGL2", "CSTA", "CYBB", "MTRNR2L12", "CTSS"),
  c("STMN1", "MKI67", "TYMS", "PCNA", "TUBA1B", "TUBB", "DUT", "HMBG2", "HMGN2", "HIST1H4C"),
  c("FCER1A", "HLA-DRB5", "HLA-BRB1"),
  c("JCHAIN", "IGHG1", "IGHG3", "MYB1", "IGHA1", "IGLC7", "IGLC2", "IGHM", "IGLC3", "IGKC"),
  c("PPBP", "PF4", "NRGN", "GNG11", "TUBB1", "CAVIN2", "HBG2", "MYL9", "GP9", "CLU"),
  c("HBA2", "HBA1", "HBD", "ALAS2", "HBM", "AHSP", "CA1", "SNCA", "HBB"),
  c("SOX4", "PRSS57", "SPINK2", "AVP", "EGFL7", "CDK6", "ANKRD28", "SERPINB1", "AREG"),
  c("ITM2C", "PLD4", "SERPINF1", "TCF4", "LILRA4", "IRF8", "IRF7", "PTGDS", "PPP1R14B"))

names(marker.list) <- c("CD4T.conventional1", "CD4T.conventional2", "Classical.Mono", "CD8T.EffectorMemory", "NK", "InflamMono", "B", "CD8T.N.NM", "Neutrophil",
                        "Non.Classical", "Immature.Neutro", "MonoZFP36L2", "Prol.T", "MDCs", "Plasmablasts", "Mega", "Erythroblast", "Undefined", "PDCs")

# Plotting the markers by individual cell types
# The markers are from our Covid-19 paper
for(i in 1:length(marker.list)){
  print(paste0("Plotting dot plot of ", names(marker.list[i]), " marker genes"))
  
  pdf(paste0(OUTPUT_DIR, "CellType/DotPlot_Covid19_", names(marker.list)[i], ".pdf"), width = 10, height = 8)
  print(DotPlot(data,features=marker.list[[i]],cols="RdBu")+coord_flip() +
          theme(legend.position = "top",legend.title = element_blank())+xlab("")+ylab("") +
          ggtitle(label = paste0("PBMC Marker ", names(marker.list)[i])))
  dev.off()
}


# We generate the violin plots for the PBMC marker genes
violin.data <- get.violin.data1(seurat = data, genes = markers_pbmc)

pdf(paste0(OUTPUT_DIR, "CellType/Violin_CellTypeMarkersPBMC.pdf"), width = 24, height = 9)
print(ggplot(violin.data, aes(factor(ident), value)) +
        geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) + 
        ylab("") + xlab("") + coord_flip() +
        facet_wrap(~ gene,scales = "free_x", ncol = length(unique(violin.data$gene))) +
        theme(strip.text.x = element_text(size=30, angle=50),
              axis.text.y = element_text(size = 30),
              strip.background = element_blank(),
              panel.spacing.x = unit(c(-0.2), "lines"),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
        scale_x_discrete(limits = rev(levels(violin.data$ident)), position = "top"))
dev.off()

# We generate the violin plot for the neutrophile marker genes
violin.data <- get.violin.data1(seurat = data, genes = markers_neutrophil)

pdf(paste0(OUTPUT_DIR, "CellType/Violin_CellTypeMarkersNeutrophile.pdf"), width = 20, height = 9)
print(ggplot(violin.data, aes(factor(ident), value)) +
        geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) +
        ylab("") + xlab("") + coord_flip() +
        facet_wrap(~ gene,scales = "free_x", ncol = length(unique(violin.data$gene))) +
        theme(strip.text.x = element_text(size=30, angle=50),
              axis.text.y = element_text(size = 30),
              strip.background = element_blank(),
              panel.spacing.x = unit(c(-0.2), "lines"),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
        scale_x_discrete(limits = rev(levels(violin.data$ident)), position = "top"))
dev.off()

# We generate the violin plot for the bmmc marker genes
violin.data <- get.violin.data1(seurat = data, genes = markers_bmmc)

pdf(paste0(OUTPUT_DIR, "CellType/Violin_CellTypeMarkersBMMC.pdf"), width = 20, height = 9)
print(ggplot(violin.data, aes(factor(ident), value)) +
        geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) +
        ylab("") + xlab("") + coord_flip() +
        facet_wrap(~ gene,scales = "free_x", ncol = length(unique(violin.data$gene))) +
        theme(strip.text.x = element_text(size=30, angle=50),
              axis.text.y = element_text(size = 30),
              strip.background = element_blank(),
              panel.spacing.x = unit(c(-0.2), "lines"),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
        scale_x_discrete(limits = rev(levels(violin.data$ident)), position = "top"))
dev.off()

# We generate the violin plot for the bone marrow marker genes
violin.data <- get.violin.data1(seurat = data, genes = markers_mnc)

pdf(paste0(OUTPUT_DIR, "CellType/Violin_CellTypeMarkersMNC.pdf"), width = 20, height = 9)
print(ggplot(violin.data, aes(factor(ident), value)) +
        geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) +
        ylab("") + xlab("") + coord_flip() +
        facet_wrap(~ gene,scales = "free_x", ncol = length(unique(violin.data$gene))) +
        theme(strip.text.x = element_text(size=30, angle=50),
              axis.text.y = element_text(size = 30),
              strip.background = element_blank(),
              panel.spacing.x = unit(c(-0.2), "lines"),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
        scale_x_discrete(limits = rev(levels(violin.data$ident)), position = "top"))
dev.off()

# We generate Dot plots of the marker genes
pdf(paste0(OUTPUT_DIR, "CellType/DotPlot_PBMC_markers.pdf"), width = 12, height = 9)
print(DotPlot(data,features=markers_pbmc,cols="RdBu")+coord_flip() +
        theme(legend.position = "top",legend.title = element_blank())+xlab("")+ylab("") +
        ggtitle(label = "PBMC Marker Genes"))
dev.off()

pdf(paste0(OUTPUT_DIR, "CellType/DotPlot_Neutrophil_markers.pdf"), width = 12, height = 9)
print(DotPlot(data,features=markers_neutrophil,cols="RdBu")+coord_flip() +
        theme(legend.position = "top",legend.title = element_blank())+xlab("")+ylab("") +
        ggtitle(label = "Neutrophil Marker Genes"))
dev.off()

pdf(paste0(OUTPUT_DIR, "CellType/DotPlot_BMMC_markers.pdf"), width = 14, height = 9)
print(DotPlot(data,features=markers_bmmc,cols="RdBu")+coord_flip() +
        theme(legend.position = "top",legend.title = element_blank())+xlab("")+ylab("") +
        ggtitle(label = "BMMC Marker Genes"))
dev.off()
#########################################################################################################
# cell types
cell.type <- data@meta.data$seurat_clusters
levels(cell.type) <- as.factor(c("Naïve CD4+ T","CD8+ T","CD14+ Monocyte","NK","Naïve CD4+ T","B","CD14+ Monocyte","Naïve CD4+ T","CD8+ T","CD16+ Monocyte",
                                 "Erythroblast","Prol. T","Prol. T","Mem. CD4+ T","Megakaryocyte","CD34+ CD38+ Progenitors","MDC","CD14+ Monocyte","PDC","Plasmablast"))


data@meta.data$BroadCelltype <- cell.type
data$celltype <- as.factor(paste0(data$BroadCelltype, " - ",data$Tissue))


pdf(paste0(OUTPUT_DIR, "CellType/UMAP_CellType.pdf"), width = 18, height = 10)
print(DimPlot(data, group.by = "celltype", pt.size = 0.5, reduction = "umap", label = T) +
        ggtitle(label = "UMAP of cell types") +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()

pdf(paste0(OUTPUT_DIR, "CellType/UMAP_CellType2.pdf"), width = 18, height = 10)
print(DimPlot(data, group.by = "celltype", pt.size = 0.5, reduction = "umap", label = T) +
        ggtitle(label = "UMAP of cell types") +
        NoLegend() +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()


pdf(paste0(OUTPUT_DIR, "CellType/UMAP_BroadCellType.pdf"), width = 18, height = 10)
print(DimPlot(data, group.by = "BroadCelltype", pt.size = 0.5, reduction = "umap", label = T) +
        ggtitle(label = "UMAP of cell types") +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()

pdf(paste0(OUTPUT_DIR, "CellType/UMAP_BroadCellType2.pdf"), width = 18, height = 10)
print(DimPlot(data, group.by = "BroadCelltype", pt.size = 0.5, reduction = "umap", label = T) +
        ggtitle(label = "UMAP of cell types") +
        NoLegend() +
        theme(text = element_text(size = 30),
              axis.text = element_text(size = 20)))
dev.off()


#########################################################################################################

UMAP_condition <- DimPlot(data, group.by = "condition", pt.size = 1, reduction = "umap", label = F) +
  ggtitle(label = "UMAP by condition") +
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 20))

pdf(paste0(OUTPUT_DIR, "CellType/Plot_Condition.pdf"), width = 12, height = 7)
print(UMAP_condition)
dev.off()

rm(marker.list, UMAP_condition, violin.data, cell.type, Female.ylinked, i, Male.xic, markers_bmmc, markers_covid19, 
   markers_mnc, markers_neutrophil, markers_pbmc, s, xic, XIC, ylinked, YLINKED_GENES)

# We save the annotated object
save.image(paste0("/Users/kuijpers/Documents/PhD/SepsisThyroid/TC/plots/Both/data/Thyroid.RData"))

