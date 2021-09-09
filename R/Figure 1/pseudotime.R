######. install.library monocle3  in R
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install()

# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor'))

# install.packages("devtools")
# devtools::install_github('cole-trapnell-lab/leidenbase')
# devtools::install_github('cole-trapnell-lab/monocle3')

#### if meet bug, try this: 
# install dependency in Mac Terminal (not need for linux)
# 1) xcode-select --install
# 2) Download new gfortran binaries for your operating system and decompress the folder (eg: gunzip gfortran-8.3-bin.tar.gz). Then do, sudo tar -xvf gfortran-8.3-bin.tar -C /. This installs everything in /usr/local/gfortran.
#######


# here data is our Seurat object with normalization and celltypes already assigned
load(paste0("/Users/kuijpers/Documents/PhD/SepsisThyroid/TC/plots/Both/data/Thyroid.RData"))

library(lattice)
library(cowplot)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(monocle3)
library(Seurat)

get_earliest_principal_node <- function(cds){
  cell_ids <- dim(cds)[2]
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

data$pool <- as.character(data$pool)
data2 <- data

for (ct in unique(data@meta.data$celltype)) {
  
  TC <- subset(data, celltype == ct)
  
  # We get the data and meta data from our Seurat Objects
  data2 <- as(as.matrix(TC@assays$RNA@data), 'sparseMatrix')
  pd <- TC@meta.data
  fData <- data.frame(gene_short_name = row.names(data2), row.names = row.names(data2))
  cds <- new_cell_data_set(data2, cell_metadata = pd, gene_metadata = fData)
  
  # The samples have to be encoded as a character otherwise, they we cannot colour the cells by sample
  cds@colData$sample <- as.character(cds@colData$sample)
  # Because we have normalized before in Seurat, we do not need to normalise here again.
  # Since we used log-normalisation, we can also subset the data, without destroying the normalisation.
  cds <- preprocess_cds(cds, num_dim = 30, norm_method = "none") 
  cds <- align_cds(cds, alignment_group = NULL, residual_model_formula_str = "~ percent.mt") #here add all unwanted co-effects
  cds <- reduce_dimension(cds)
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  
  pdf(paste0(OUTPUT_DIR,"monocle3/",ct,"_Trajectory_Condition.pdf"), width = 6, height = 5)
  print(plot_cells(cds, color_cells_by = "condition", graph_label_size = 5, 
                   label_leaves = F, label_branch_points = F, cell_size = 0.75, trajectory_graph_segment_size = 1,
                   label_cell_groups = F, label_groups_by_cluster = F) + 
          theme(legend.position = "right"))
  dev.off()
  
  pdf(paste0(OUTPUT_DIR,"monocle3/",ct,"_Trajectory_Sample.pdf"), width = 6, height = 5)
  print(plot_cells(cds, color_cells_by = "sample", graph_label_size = 5, 
                   label_leaves = F, label_branch_points = F, cell_size = 0.75, trajectory_graph_segment_size = 1,
                   label_cell_groups = F, label_groups_by_cluster = F) + 
          theme(legend.position = "right"))
  dev.off()
  
  pdf(paste0(OUTPUT_DIR,"monocle3/",ct,"_Trajectory_batch.pdf"), width = 6, height = 5)
  print(plot_cells(cds, color_cells_by = "pool", graph_label_size = 5, 
                   label_leaves = F, label_branch_points = F, cell_size = 0.75, trajectory_graph_segment_size = 1,
                   label_cell_groups = F, label_groups_by_cluster = F) + 
          theme(legend.position = "right"))
  dev.off()
  
  cds.ori <- cds
  
  # To order the cells according to their pseudo time, we have to determine a start point.
  # Since we have no data indicating which cells are older, we will determine the start point programmatically.
  # a helper function to identify the root principal points
  # Each cell has a closest root. We choose the root, that is the closest root for the most cells. 
  
  cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
  trace1 <- plot_cells(cds, label_cell_groups = F, color_cells_by = "pseudotime", label_leaves = F, 
                       label_branch_points = F, graph_label_size = 5, cell_size=0.5, trajectory_graph_segment_size = 2)
  
  pdf(paste0(OUTPUT_DIR,'monocle3/',ct,'_PseudoTime.pdf'), width = 6, height = 5)
  print(trace1)
  dev.off()
}

############################################################################################################################################
# Here we plot the cells using a custom setup. We basically use parts of the plot_cells funtion and 
# supply our own data.frame.
# We get the UMAP values per cell and the meta data we colour them by. 
trajectory_graph_segment_size = 1
trajectory_graph_color = "gray29"
color_cells_by = "condition"
cell_size = 0.2

custom.plot.data <- data.frame(reducedDims(cds)[["UMAP"]], 
                               condition = factor(cds@colData$condition, levels = c("Healthy", "Goiter", "TC")),
                               celltype = cds@colData$celltype,
                               pseudotime = pseudotime(cds, reduction_method = "UMAP")
)
colnames(custom.plot.data)[1:2] <- c("UMAP1", "UMAP2")

# We also want to show the trajectory and get the relevant information out of the cds object
ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>% 
  as.data.frame() %>% 
  dplyr::select(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>% 
  dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))

dp_mst <- cds@principal_graph[["UMAP"]]
edge_df <- dp_mst %>% 
  igraph::as_data_frame() %>% 
  dplyr::select(source = "from", target = "to") %>% 
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(source = "sample_name", 
                                   source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                   source_prin_graph_dim_2 = "prin_graph_dim_2"), 
                   by = "source") %>% 
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(target = "sample_name", 
                                   target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                   target_prin_graph_dim_2 = "prin_graph_dim_2"), 
                   by = "target")


pdf(paste0(OUTPUT_DIR,"monocle3/Condition_Trajectory_Condition_custom.pdf"), width = 6, height = 5)
ggplot(custom.plot.data, aes(x = UMAP1, y = UMAP2, colour = condition)) +
  geom_point(size = cell_size) +
  theme_classic() +
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = trajectory_graph_segment_size, 
               color = I(trajectory_graph_color), linetype = "solid", na.rm = TRUE, data = edge_df) +
  guides(color = guide_legend(title = color_cells_by, 
                              override.aes = list(size = 4)))
dev.off()
############################################################################################################################################


# We plot the change of the top DEGs in the trajectory
degs <- read.csv(paste0(OUTPUT_DIR,"/DEG/DEG_bulk_verysig.csv"), stringsAsFactors = FALSE)
degs.up <- subset(degs, avg_logFC > 0)
degs.down <- subset(degs, avg_logFC < 0)

# We now test which genes are changing expression along the trajectory
cds.pr.test <- graph_test(cds, neighbor_graph = "principal_graph", cores = 7)
cds.pr.test <- cds.pr.test[order(cds.pr.test$p_value),]
cds.pr.test$p_value <- p.adjust(cds.pr.test$p_value, method = "BH")
cds.pr.test.sig <- subset(cds.pr.test, p_value < 0.01)

# We read in the list of GWAS sepsis genes 
#gwas.genes <- read.csv("/Users/grasshoff/Documents/SepsisThyroidAnalysis/reference/gwas.genes.csv")

#cds.pr.test.sig.sepsis <- subset(cds.pr.test.sig, gene_short_name %in% gwas.genes$genes)

#pdf("/Users/grasshoff/Documents/SepsisThyroidAnalysis/output/scRNA/monocle3/Condition_Trajectory_Sepsis.pdf", width = 12, height = 10)
#plot_cells(cds,
#           genes=cds.pr.test.sig.sepsis$gene_short_name[1:dim(cds.pr.test.sig.sepsis)[1]],
#           label_cell_groups=TRUE,
#           show_trajectory_graph=TRUE, 
#           graph_label_size = 1, 
#           cell_size=1, 
#           trajectory_graph_segment_size = 2,
#           label_leaves = F, label_branch_points = F
#           )
#dev.off()


# Test if any of the DEGs are significantly changing over the course of the trajectory
cds.pr.test.sig.degs <- subset(cds.pr.test.sig, gene_short_name %in% degs$gene)
intersected <- sum(degs$gene %in% cds.pr.test.sig.degs$gene)
title.text <- paste0("All cells: ", dim(degs)[1], " DEGs\n", intersected, " DEGs show a significant change along the trajectory")

pdf(paste0(OUTPUT_DIR,"monocle3/Condition_Trajectory_DEGsTop9.pdf"), width = 9, height = 7.5)
plot_cells(cds,
           genes=cds.pr.test.sig.degs$gene_short_name[1:9],
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE, 
           graph_label_size = 5, 
           cell_size=1.5, 
           trajectory_graph_segment_size = 2,
           label_leaves = F, label_branch_points = F
) +
  ggtitle(label = title.text)
dev.off()

# We write the the results 
cds.pr.test.sig.degs <- data.frame(gene = rownames(cds.pr.test.sig.degs), cds.pr.test.sig.degs)
write.csv(cds.pr.test.sig.degs, paste0(OUTPUT_DIR,"monocle3/Condition_Trajectory_DEGs.csv"), quote = FALSE, row.names = FALSE)

# Test which of the up regulated genes change along the trajectory
cds.pr.test.sig.degs.up <- subset(cds.pr.test.sig, gene_short_name %in% degs.up$gene)
intersected <- sum(degs.up$gene %in% cds.pr.test.sig.degs.up$gene)
title.text <- paste0("All cells: ", dim(degs.up)[1], " up regulated DEGs\n", intersected, " DEGs show a significant change along the trajectory")

pdf(paste0(OUTPUT_DIR,"monocle3/Condition_Trajectory_DEGs_Up_Top9.pdf"), width = 9, height = 7.5)
plot_cells(cds,
           genes=cds.pr.test.sig.degs.up$gene_short_name[1:9],
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE, 
           graph_label_size = 5, 
           cell_size=1.5, 
           trajectory_graph_segment_size = 2,
           label_leaves = F, label_branch_points = F
) +
  ggtitle(label = title.text)
dev.off()

# Test which of the down regulated genes change along the trajectory
cds.pr.test.sig.degs.down <- subset(cds.pr.test.sig, gene_short_name %in% degs.down$gene)
intersected <- sum(degs.down$gene %in% cds.pr.test.sig.degs.down$gene)
title.text <- paste0("All cells: ", dim(degs.down)[1], " down regulated DEGs\n", intersected, " DEGs show a significant change along the trajectory")

pdf(paste0(OUTPUT_DIR,"monocle3/Condition_Trajectory_DEGs_Down_Top9.pdf"), width = 9, height = 7.5)
plot_cells(cds,
           genes=cds.pr.test.sig.degs.down$gene_short_name[1:9],
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE, 
           graph_label_size = 5, 
           cell_size=1.5, 
           trajectory_graph_segment_size = 2,
           label_leaves = F, label_branch_points = F
) +
  ggtitle(label = title.text)
dev.off()