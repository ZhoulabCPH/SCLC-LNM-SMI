
# Figure 2A
library(Seurat)
library(tidyverse)
library(ggplot2)
library(harmony)
library(ROGUE)
library(ggpubr)
malignant$type <- paste(malignant$Patient_metastasis,
                        malignant$Tissue_metastasis)
malignant_PT <- subset(malignant,subset = type == "Primary Primary")
malignant_PT_LNM <- subset(malignant,subset = type == "Metastasis Primary")
malignant_LNMT <- subset(malignant,subset = type == "Metastasis LNM")
# number double check
length(unique(substr(malignant_PT$TMA,1,7)))
length(unique(substr(malignant_PT_LNM$TMA,1,7)))
length(unique(substr(malignant_LNMT$TMA,1,7)))
##  total patient
dim_plot <- DimPlot(malignant, reduction = "umap", 
                    cols = c("#AF8260","#F8B195",
                             "#F39C12","#D2691E","#FF6347","#C63C51",
                             "#E25A90","#EEC60A","#BB8493","#FF8D00",
                             "#B8860B","#F4ACB7","#8D493A"),
                    group.by = "seurat_clusters",raster=T) 
LabelClusters(dim_plot, id = "seurat_clusters", repel = TRUE)
##  PT
DimPlot(malignant_PT, reduction = "umap",  
        cols = c("#AF8260","#F8B195",
                 "#F39C12","#D2691E","#FF6347","#C63C51",
                 "#E25A90","#EEC60A","#BB8493","#FF8D00",
                 "#B8860B","#F4ACB7","#8D493A"),
        group.by = "seurat_clusters",raster=T) 
##  PT-LNM
DimPlot(malignant_PT_LNM, reduction = "umap",  
        cols = c("#AF8260","#F8B195",
                 "#F39C12","#D2691E","#FF6347","#C63C51",
                 "#E25A90","#EEC60A","#BB8493","#FF8D00",
                 "#B8860B","#F4ACB7","#8D493A"),
        group.by = "seurat_clusters",raster=T) 
##  LNMT
DimPlot(malignant_LNMT, reduction = "umap", 
        cols = c("#AF8260","#F8B195",
                 "#F39C12","#D2691E","#FF6347","#C63C51",
                 "#E25A90","#EEC60A","#BB8493","#FF8D00",
                 "#B8860B","#F4ACB7","#8D493A"),
        group.by = "seurat_clusters",raster=T) 


# Figure 2B
library(dplyr)
options(stringsAsFactors=FALSE)
library(reticulate)
divMatrix <- function(mat, divisor) {
  return(mat / divisor)
}
# +++, Ro/e > 3; ++, 1 <Ro/e ≤ 3; +, 0.2 ≤ Ro/e ≤ 1; +/−, 0 < Ro/e < 0.2; −, Ro/e = 0.
ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}
ROIE_with_pvalue <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab and p-values
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix with p-values
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  p_values <- matrix(NA, nrow = nrow(crosstab), ncol = ncol(crosstab))
  for (i in 1:nrow(crosstab)) {
    for (j in 1:ncol(crosstab)) {
      observed <- crosstab[i, j]
      expected <- (rowSums(crosstab)[i] * colSums(crosstab)[j]) / allsum
      test <- fisher.test(rbind(c(observed, rowSums(crosstab)[i] - observed),
                                c(colSums(crosstab)[j] - observed, allsum - rowSums(crosstab)[i] - colSums(crosstab)[j] + observed)))
      p_values[i, j] <- test$p.value
    }
  }
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  row.names(p_values) <- row.names(crosstab)
  colnames(p_values) <- colnames(crosstab)
  return(list(RoE = roie, p_values = p_values))
}
library(pheatmap)
malignant$type <- paste(malignant$Patient_metastasis,
                        malignant$Tissue_metastasis)
malignant_PT <- subset(malignant,subset = type == "Primary Primary")
malignant_PT_LNM <- subset(malignant,subset = type == "Metastasis Primary")
malignant_LNMT <- subset(malignant,subset = type == "Metastasis LNM")
malignant_PT_id <- table(malignant_PT$ID,
                         malignant_PT$seurat_clusters)
malignant_PT_LNM_id <-  table(malignant_PT_LNM$ID,
                              malignant_PT_LNM$seurat_clusters)
malignant_LNMT_id <-  table(malignant_LNMT$ID,
                            malignant_LNMT$seurat_clusters)
malignant_mat <- data.frame(PT = colSums(malignant_PT_id),
                            PT_LNM = colSums(malignant_PT_LNM_id),
                            LNMT = colSums(malignant_LNMT_id))
malignant_mat_roie1 <- ROIE(malignant_mat)

malignant_mat_roie <- ROIE_with_pvalue(malignant_mat)
malignant_mat_roie_P <- malignant_mat_roie$p_values
malignant_mat_roie <- malignant_mat_roie$RoE
library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(c(0,1),c("#000000","#FFD700"))
malignant_mat_roie <- as.matrix(malignant_mat_roie)
Heatmap(malignant_mat_roie, name ="malignant_mat_roie",col = col_fun,
        cluster_columns = F,cluster_rows = F,
        width = unit(12, "mm") * ncol(malignant_mat_roie),
        height = unit(10, "mm") * nrow(malignant_mat_roie),
        cell_fun=function(j,i,x,y,width, height,fill)
        {grid.text(sprintf("%.2f",malignant_mat_roie[i,j]),
                   x,y,gp= gpar(fontsize = 10))})
col_fun <- colorRamp2(c(0,1),c("#FF6F61","#A9B7C6"))
Heatmap(malignant_mat_roie_P, name ="malignant_mat_roie_P",col = col_fun,
        cluster_columns = F,cluster_rows = F,
        width = unit(12, "mm") * ncol(malignant_mat_roie_P),
        height = unit(10, "mm") * nrow(malignant_mat_roie_P),
        cell_fun=function(j,i,x,y,width, height,fill)
        {grid.text(sprintf("%.3f",malignant_mat_roie_P[i,j]),
                   x,y,gp= gpar(fontsize = 10))})

# Figure 2C
library(Seurat)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(grid)
library(gghalves)
library(ggpubr)
malignant_marker <- malignant_marker[which(malignant_marker$avg_log2FC > 1.1),]
malignant_marker$cluster <- factor(malignant_marker$cluster,
                                   levels = c("Malignant C1","Malignant C2",
                                              "Malignant C3","Malignant C4",
                                              "Malignant C5","Malignant C6",
                                              "Malignant C7","Malignant C8",
                                              "Malignant C9","Malignant C10",
                                              "Malignant C11","Malignant C12",
                                              "Malignant C13"))
markers <- c(malignant_marker$gene)
uni_markers <- names(table(markers)[which(table(markers) < 2)])
malignant_marker <- malignant_marker[match(uni_markers,malignant_marker$gene),]
malignant_marker <- malignant_marker[order(malignant_marker$cluster),]
malignant_data <- AverageExpression(malignant,
                                    assays = "RNA", 
                                    return.seurat = TRUE,
                                    group.by = "clusters",
                                    slot = "data")
cluster.averages <- ScaleData(malignant_data)
mat <- GetAssayData(cluster.averages, slot = "scale.data")
mat <- as.matrix(mat[malignant_marker$gene, ])
pheatmap(t(mat),fontsize=6,
         color  = colorRampPalette(c("#7F3C8D","#F2F2F2","#FFCB05"))(100),
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
# Malignant C5
ggplot(bp_c5_result, aes(x=Description, y=-log10(p.adjust),
                         fill = "#FF6347")) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("-log10(FDR)") +
  theme(axis.text.x=element_text(hjust = 1,colour="black",
                                 size=10,angle = 90))
# Malignant C6
ggplot(bp_c6_result, aes(x=Description, y=-log10(p.adjust),
                         fill = "#C63C51")) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("-log10(FDR)") +
  theme(axis.text.x=element_text(hjust = 1,colour="black",
                                 size=10,angle = 90))
# Malignant C9
ggplot(bp_c9_result, aes(x=Description, y=-log10(p.adjust),
                         fill = "#BB8493")) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("-log10(FDR)") +
  theme(axis.text.x=element_text(hjust = 1,colour="black",
                                 size=10,angle = 90))

























