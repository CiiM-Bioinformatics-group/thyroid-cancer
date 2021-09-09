# data <- readRDS('../data.RDS')
library(Seurat)
library(ggplot2)

rm(list = ls())

data <- subset(data, celltype != 'B - PBMC')
data@meta.data[which(data@meta.data$celltype == 'CLP'), 'celltype'] <- 'Lymphoid progenitors'
data@meta.data[which(data@meta.data$celltype == 'MPP'), 'celltype'] <- 'Myeloid/multipotent progenitors'
data@meta.data[which(data@meta.data$celltype == 'unclassified'), 'celltype'] <- 'Early B cell progenitors'
data@meta.data[which(data@meta.data$celltype == 'Pre-B Cell'), 'celltype'] <- 'Late B cell progenitors'

data@meta.data[which(data$condition == 'MNS'), ]$condition <- 'MNG'


TestingCellProportions <- function(data, identity.to.compare, condition.to.compare, output_path, type, celltypes, PlotSigIndicators = TRUE){
  
  # Libraries
  library(Seurat)
  library(ggplot2)
  library(ggsignif)
  library(reshape2)
  if(!require("DirichletReg")) install.packages("DirichletReg")
  
  # Reshaping the data for the Dirichlet-Test
  # Descriptive Statistics
  # Mean, Variance and Coefficient of Variance per Cell Type
  DescriptionCelltype <- function(SeuratObject, identity.to.compare, condition.to.compare, celltypes){
    
    # Result object
    result <- list()
    
    # Subsetting to the relevant subgroup
    data <- FetchData(SeuratObject, vars = identity.to.compare)
    data.subgroup <- SeuratObject[, which(x = data == condition.to.compare)]
    
    # Calculating the percentages per cell type per group
    result$CellsPerType <- table(as.character(data.subgroup$celltype))
    result$PercentagesPerType <- result$CellsPerType/dim(data.subgroup)[2]
    
    # Calculating the percentages per cell type per patient
    patients <- unique(data.subgroup$sample)
    patients <- patients[order(patients)]
    print(patients)
    CellsPerTypePerPatient <- c()
    PercentagesPerTypePerPatient <- c()
    values.for.MeanVariance <- c()
    for(person in patients){
      # Not all patients have all cell types at all time points.
      # With this dummy, we can use 0 for missing values
      result.dummy <- rep(0,length(celltypes))
      names(result.dummy) <- celltypes
      result.dummy <- result.dummy[sort(names(result.dummy))]
      
      # We subset to a sample ID (one person at one time point)
      data.person <- FetchData(data.subgroup, vars = "sample")
      data.person <- data.subgroup[, which(x = data.person == person)]
      
      # How many cells does a patient have of each type?
      result.preliminary <- table(as.character(data.person[[type]][,1]))
      result.preliminary <- result.preliminary[sort(names(result.preliminary))]
      
      result.dummy[names(result.dummy) %in% names(result.preliminary)] <- result.preliminary
      CellsPerTypePerPatient <- rbind(CellsPerTypePerPatient, result.dummy)
      rownames(CellsPerTypePerPatient)[dim(CellsPerTypePerPatient)[1]] <- paste0("Patient_", person)
      
      
      # Percentages
      result.dummy <- result.dummy/sum(result.dummy)
      PercentagesPerTypePerPatient <- rbind(PercentagesPerTypePerPatient, result.dummy)
      rownames(PercentagesPerTypePerPatient)[dim(PercentagesPerTypePerPatient)[1]] <- paste0("Patient_", person)
    }
    
    result$CellsPerPatient <- CellsPerTypePerPatient
    result$PercentagesPerPatient <- PercentagesPerTypePerPatient
    
    # Calculating the Mean, Variance and Coefficient of Variance per cell type
    result$MeanPerCelltype <- colMeans(PercentagesPerTypePerPatient)
    result$VarsPerCelltype <- apply(PercentagesPerTypePerPatient, 2, var)
    result$CoeffOfVariance <- result$VarsPerCelltype / result$MeanPerCelltype
    result$Condition <- condition.to.compare
    return(result)
  }
  
  CellTypePerCondition <- lapply(condition.to.compare, DescriptionCelltype,
                                 SeuratObject = data,
                                 identity.to.compare = identity.to.compare,
                                 celltypes = celltypes)
  
  print(CellTypePerCondition)
  
  # Boxplot to illustrate the variance
  data.boxplot.combined <- c()
  for(conditions in 1:length(condition.to.compare)){
    data.boxplot <- data.frame(patient = rownames(CellTypePerCondition[[conditions]]$PercentagesPerPatient), CellTypePerCondition[[conditions]]$PercentagesPerPatient, stringsAsFactors = FALSE)
    data.boxplot <- melt(data.boxplot, id.vars = "patient")
    
    # pdf(paste0(output_path, "Condition_", condition.to.compare[conditions], "_CellTypes.pdf"))
    # print(
    #   ggplot(data.boxplot, aes(x = variable, y = value)) +
    #     geom_boxplot() +
    #     theme(axis.text.x = element_text(angle = 45, hjust = 1)))
    # dev.off()
    #
    # Combining the data to plot all conditions in one boxplot
    data.boxplot <- data.frame(data.boxplot, condition = condition.to.compare[conditions])
    data.boxplot.combined <- rbind(data.boxplot.combined, data.boxplot)
  }
  
  data.boxplot.combined$variable <- gsub(pattern = "\\.", replacement = ' ', x = data.boxplot.combined$variable)
  pdf(paste0(output_path, "AllConditions_CellTypes.pdf"), width = 10, height = 4)
  print(
    ggplot(data = data.boxplot.combined, aes(x = condition, y = value, fill = condition))+
      geom_boxplot(outlier.shape = NA)+
      facet_grid(.~variable)+
      geom_jitter(shape=16, width = 0.2)+
      theme_classic() +
      ylab("Cell proportions" )+
      scale_fill_manual(values=c("#00AFBB",
                                 "#E7B800",
                                 "#FC4E07"))+
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  )
  dev.off()
  
  # Performing the Dirichlet-Tests for each Cell Type
  # We loop over all conditions by prepending them with a LETTER
  
  results.all <- c() # We rbind our results to this variable
  for(i in 1:length(as.character(unique(data[[identity.to.compare]])[,1]))){
    # What are the LETTERS with which we prepend the conditions
    prepend <- LETTERS[unique(c(i:length(as.character(unique(data[[identity.to.compare]])[,1])), 1:i))]
    print(cat("What we prepend with: ", prepend))
    
    data.test <- data.frame(celltype  = as.character(data.boxplot.combined$variable),
                            sample    = as.character(data.boxplot.combined$patient),
                            value     = data.boxplot.combined$value,
                            other     = 1 - data.boxplot.combined$value,
                            condition = as.character(data.boxplot.combined$condition),
                            stringsAsFactors = FALSE)
    celltypes = unique(data.test$celltype)
    
    conditions <- unique(data.test$condition)
    conditions <- conditions[order(conditions)]
    
    for(prepending in 1:length(prepend)){
      data.test[data.test$condition == conditions[prepending], 5] <- paste0(prepend[prepending],conditions[prepending])
    }
    
    # Vector of all states in the identity to compare
    condition.vector <- unique(data.test$condition)
    condition.vector <- condition.vector[order(condition.vector)]
    comparisons <- paste0("condition", condition.vector)
    dat.test <- NULL
    for(type in celltypes){
      try({
        print(type)
        data.type <- data.test[data.test$celltype == type,]
        data.type <- rbind(subset(data.type, condition == condition.vector[1]),
                           subset(data.type, condition == condition.vector[2]),
                           subset(data.type, condition == condition.vector[3]),
                           subset(data.type, condition == condition.vector[4]))
        data.type$Smp <- DR_data(data.type[,3:4])
        res <- DirichReg(Smp~condition, data.type, model = "alternative", base = 1, verbosity = 0, control = list(iterlim = 5000))
        x <- summary(res)
        
        for(comparison in comparisons){
          if(comparison %in% rownames(x$coef.mat)){
            dat.test <- rbind(dat.test, c(type, x$coef.mat[comparison, 3:4], comparison))
            # Renaming the comparisons
            dat.test[dat.test[,4] == comparison,4] <- paste0(substr(condition.vector[1], start = 2, stop = nchar(condition.vector[1])), "Vs", substr(comparison, start = 11, stop = nchar(comparison)))
          }
        }
      })
    }
    
    dat.test <- as.data.frame(dat.test, stringsAsFactors = FALSE)
    colnames(dat.test) <- c("cell_type","Z","Pr","comparison")
    
    # Saving all results
    results.all <- rbind(results.all, dat.test)
  }
  # Renaming all the comparisons and removing duplicated ones
  old.comparisons <- results.all$comparison
  old.comparisons <- strsplit(old.comparisons, "Vs")
  old.comparisons <- lapply(old.comparisons, sort)
  new.comparisons <- lapply(old.comparisons, paste0, collapse = "Vs")
  new.comparisons <- unlist(new.comparisons)
  results.all$comparison <- paste0(results.all$cell_type, new.comparisons)
  results.all <- results.all[!duplicated(results.all$comparison),]
  
  # Bonferroni correcting the p-values
  results.all$sig <- "ns"
  results.all$PrAdjust <- p.adjust(results.all$Pr, method = "bonferroni")
  results.all$sig[results.all$PrAdjust < 0.05]="*"
  results.all$sig[results.all$PrAdjust < 0.01]="**"
  results.all$sig[results.all$PrAdjust < 0.001]="***"
  
  # Writing all results
  write.csv(results.all, paste0(output_path, "DirichletTestAdjusted.csv"), row.names = FALSE)
  # Saving only results that significant
  results.all.sig <- subset(results.all, sig != "ns")
  write.csv(results.all.sig, paste0(output_path, "DirichletTestAdjustedSig.csv"), row.names = FALSE)
}

identity.to.compare <- "condition"
condition.to.compare <- as.character(levels(as.factor(data@meta.data$condition)))
type <- "celltype"
celltypes <- unique(as.character(data@meta.data$celltype))
OUTPUT_DIR <- '/vol/projects/mzoodsma/thyroid_cancer/output/cell_props/'

output_path <- paste0(OUTPUT_DIR, "broad/")
TestingCellProportions(data,
                       identity.to.compare,
                       condition.to.compare,
                       output_path,
                       type,
                       celltypes,
                       PlotSigIndicators = TRUE)