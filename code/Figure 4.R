###   PT，PT-LNM 和 LNMT 分开算 ###
# y_local_px: Same as “x_local_px” but for the y dimension. 
# Due to differences in Cartesian and imaging plotting frames of reference, 
# the local y-axis coordinates should be inverted to align the 
# transcript positions to images.
# To convert to microns multiply the pixel value by 0.12 um per pixel
library(dbscan)
library(foreach)
library(doParallel)
PT <- readRDS("PT_harmony8_afterQC_Doubcell_result.RDS")
PT_LNM <- readRDS("PT_LNM_harmony8_afterQC_Doubcell_result.RDS")
LNMT <- readRDS("LNMT_harmony8_afterQC_Doubcell_result.RDS")
PT_TMA <- names(table(PT$TMA))
PT_LNM_TMA <- names(table(PT_LNM$TMA))
LNMT_TMA <- names(table(LNMT$TMA))
numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)
### PT
score_PT <- list()
Pvalue_PT <- list()
for (q in 1:length(PT_TMA)){
  weizhi <- which(PT$TMA == PT_TMA[q])
  spatial_mat <- data.frame(x_μm = PT$x_slide_mm[weizhi] * 1000,
                            y_μm = PT$y_slide_mm[weizhi] * 1000,
                            cell = PT$cluster[weizhi])
  coords <- as.matrix(spatial_mat[, c("x_μm", "y_μm")])
  threshold <- 6
  neighbor_list <- frNN(coords, eps = threshold)
  cell_types <- unique(spatial_mat$cell)
  observed_interactions <- matrix(0, nrow=length(cell_types), ncol=length(cell_types), 
                                  dimnames=list(cell_types, cell_types))
  cell_type_indices <- lapply(cell_types, function(type) which(spatial_mat$cell == type))
  names(cell_type_indices) <- cell_types
  for (i in seq_along(cell_types)) {
    type1 <- cell_types[i]
    indices1 <- cell_type_indices[[type1]]
    for (j in seq_along(cell_types)) {
      type2 <- cell_types[j]
      indices2 <- cell_type_indices[[type2]]
      count <- sum(sapply(indices1, function(idx) {
        neighbors <- neighbor_list$id[[idx]]
        sum(neighbors %in% indices2)
      }))
      observed_interactions[i, j] <- count
    }
  }
  num_simulations <- 1000
  random_interactions_list <- foreach(sim = 1:num_simulations, .combine = 'c', 
                                      .packages = 'dbscan') %dopar% {
                                        
                                        randomized_positions <- spatial_mat
                                        randomized_positions$x_μm <- sample(randomized_positions$x_μm)
                                        randomized_positions$y_μm <- sample(randomized_positions$y_μm)
                                        
                                        coords_rand <- as.matrix(randomized_positions[, c("x_μm", "y_μm")])
                                        neighbor_list_rand <- frNN(coords_rand, eps = threshold)
                                        
                                        sim_interactions <- matrix(0, nrow=length(cell_types), ncol=length(cell_types))
                                        for (i in seq_along(cell_types)) {
                                          type1 <- cell_types[i]
                                          indices1 <- which(randomized_positions$cell == type1)
                                          for (j in seq_along(cell_types)) {
                                            type2 <- cell_types[j]
                                            indices2 <- which(randomized_positions$cell == type2)
                                            count <- sum(sapply(indices1, function(idx) {
                                              neighbors <- neighbor_list_rand$id[[idx]]
                                              sum(neighbors %in% indices2)
                                            }))
                                            sim_interactions[i, j] <- count
                                          }
                                        }
                                        list(sim_interactions)
                                      }
  
  random_interactions <- array(unlist(random_interactions_list), 
                               dim = c(length(cell_types), length(cell_types), num_simulations))
  
  expected_interactions <- apply(random_interactions, c(1, 2), mean)
  sd_interactions <- apply(random_interactions, c(1, 2), sd)
  # interaction score
  interaction_score <- (observed_interactions - expected_interactions) / sd_interactions
  # P value
  p_values <- matrix(0, nrow=length(cell_types), ncol=length(cell_types), 
                     dimnames=list(cell_types, cell_types))
  for (i in seq_along(cell_types)) {
    for (j in seq_along(cell_types)) {
      p_values[i, j] <- sum(random_interactions[i, j, ] >= observed_interactions[i, j]) / num_simulations
    }
  }
  
  
  
  
  score_PT[[q]] <- interaction_score
  Pvalue_PT[[q]] <- p_values
}


# Figure 4B
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(metap)
all_harmony8 <- readRDS("all_harmony8_afterQC_Doubcell_result.RDS")
cell <- names(table(all_harmony8$cluster))
score_PT <- readRDS("PT_interaction_score.RDS")
score_PT_LNM <- readRDS("PT_LNM_interaction_score.RDS")
score_LNMT <- readRDS("LNMT_interaction_score.RDS")
p_PT_inter <- readRDS("PT_inter_Pvalue.RDS")
p_PT_LNM_inter <- readRDS("PT_LNM_inter_Pvalue.RDS")
p_LNMT_inter <- readRDS("LNMT_inter_Pvalue.RDS")
score_PT_full <- matrix(0, nrow=length(cell), ncol=length(cell),
                        dimnames=list(cell, cell))
score_PT_count <- matrix(0, nrow=length(cell), ncol=length(cell),
                         dimnames=list(cell, cell))
score_PT_LNM_full <- matrix(0, nrow=length(cell), ncol=length(cell),
                            dimnames=list(cell, cell))
score_PT_LNM_count <- matrix(0, nrow=length(cell), ncol=length(cell),
                             dimnames=list(cell, cell))
score_LNMT_full <- matrix(0, nrow=length(cell), ncol=length(cell),
                          dimnames=list(cell, cell))
score_LNMT_count <- matrix(0, nrow=length(cell), ncol=length(cell),
                           dimnames=list(cell, cell))
P_PT_full <- matrix(0, nrow=length(cell), ncol=length(cell),
                    dimnames=list(cell, cell))
P_PT_LNM_full <- matrix(0, nrow=length(cell), ncol=length(cell),
                        dimnames=list(cell, cell))
P_LNMT_full <- matrix(0, nrow=length(cell), ncol=length(cell),
                      dimnames=list(cell, cell))
##  merge score 
for (i in 1:length(score_PT)){
  linshi <- score_PT[[i]]
  linshi[which(linshi == "NaN")] <- NA
  score_PT_full1 <- matrix(NA, nrow=length(cell), ncol=length(cell),
                           dimnames=list(cell, cell))
  score_PT_full1[rownames(linshi), colnames(linshi)] <- as.matrix(linshi)
  score_PT_count <- score_PT_count + (!is.na(score_PT_full1))
  score_PT_full1[which(is.na(score_PT_full1) == T)] <- 0
  score_PT_full <- score_PT_full + score_PT_full1
}
score_PT_full <- score_PT_full/(score_PT_count++0.0001)
for (i in 1:length(score_PT_LNM)){
  linshi <- score_PT_LNM[[i]]
  linshi[which(linshi == "NaN")] <- NA
  score_PT_LNM_full1 <- matrix(NA, nrow=length(cell), ncol=length(cell),
                               dimnames=list(cell, cell))
  score_PT_LNM_full1[rownames(linshi), colnames(linshi)] <- as.matrix(linshi)
  score_PT_LNM_count <- score_PT_LNM_count + (!is.na(score_PT_LNM_full1))
  score_PT_LNM_full1[which(is.na(score_PT_LNM_full1) == T)] <- 0
  score_PT_LNM_full1[which(score_PT_LNM_full1 == Inf)] <- 0
  score_PT_LNM_full <- score_PT_LNM_full + score_PT_LNM_full1
}
score_PT_LNM_full <- score_PT_LNM_full/(score_PT_LNM_count+0.0001)
for (i in 1:length(score_LNMT)){
  linshi <- score_LNMT[[i]]
  linshi[which(linshi == "NaN")] <- NA
  score_LNMT_full1 <- matrix(NA, nrow=length(cell), ncol=length(cell),
                             dimnames=list(cell, cell))
  score_LNMT_full1[rownames(linshi), colnames(linshi)] <- as.matrix(linshi)
  score_LNMT_count <- score_LNMT_count + (!is.na(score_LNMT_full1))
  score_LNMT_full1[which(is.na(score_LNMT_full1) == T)] <- 0
  score_LNMT_full1[which(score_LNMT_full1 == Inf)] <- 0
  score_LNMT_full <- score_LNMT_full + score_LNMT_full1
}
score_LNMT_full <- score_LNMT_full/(score_LNMT_count+0.0001)
##   merge P  
Pvalue_PT <- list()
for (i in 1:length(score_PT)){
  linshi_score <- score_PT[[i]]
  linshi_P_inter <- p_PT_inter[[i]]
  linshi_score[which(linshi_score == "NaN")] <- 0
  linshi_P_inter[linshi_score < 0] <- NA
  linshi_P_inter[linshi_score == 0] <- NA
  P_combined <- linshi_P_inter
  Pvalue_PT[[i]] <- P_combined
}
Pvalue_PT_LNM <- list()
for (i in 1:length(score_PT_LNM)){
  linshi_score <- score_PT_LNM[[i]]
  linshi_P_inter <- p_PT_LNM_inter[[i]]
  linshi_score[which(linshi_score == "NaN")] <- 0
  linshi_P_inter[linshi_score < 0] <- NA
  linshi_P_inter[linshi_score == 0] <- NA
  P_combined <- linshi_P_inter
  Pvalue_PT_LNM[[i]] <- P_combined
}
Pvalue_LNMT <- list()
for (i in 1:length(score_LNMT)){
  linshi_score <- score_LNMT[[i]]
  linshi_P_inter <- p_LNMT_inter[[i]]
  linshi_score[which(linshi_score == "NaN")] <- 0
  linshi_P_inter[linshi_score < 0] <- NA
  linshi_P_inter[linshi_score == 0] <- NA
  P_combined <- linshi_P_inter
  Pvalue_LNMT[[i]] <- P_combined
}
p_PT <- Pvalue_PT
p_PT_LNM <- Pvalue_PT_LNM
p_LNMT <- Pvalue_LNMT
## merge P 
for (i in 1:28){
  for (j in 1:28){
    p_linshi <- c()
    for (q in 1:length(p_PT)){
      linshi <- p_PT[[q]]
      linshi[which(linshi == "NaN")] <- NA
      linshi1 <- matrix(NA, nrow=length(cell), ncol=length(cell),
                        dimnames=list(cell, cell))
      linshi1[rownames(linshi), colnames(linshi)] <- as.matrix(linshi)
      p_linshi <- c(p_linshi,
                    linshi1[rownames(P_PT_full)[i],colnames(P_PT_full)[j]])
    }
    p_linshi <- p_linshi[which(is.na(p_linshi) == F)]
    if (length(p_linshi) > 1){
      P_PT_full[i,j] <- sumlog(p = p_linshi)$p
    }
    else if (length(p_linshi) == 1){
      P_PT_full[i,j] <- p_linshi
    }
  }
}
for (i in 1:28){
  for (j in 1:28){
    p_linshi <- c()
    for (q in 1:length(p_PT_LNM)){
      linshi <- p_PT_LNM[[q]]
      linshi[which(linshi == "NaN")] <- NA
      linshi1 <- matrix(NA, nrow=length(cell), ncol=length(cell),
                        dimnames=list(cell, cell))
      linshi1[rownames(linshi), colnames(linshi)] <- as.matrix(linshi)
      p_linshi <- c(p_linshi,
                    linshi1[rownames(P_PT_LNM_full)[i],colnames(P_PT_LNM_full)[j]])
    }
    p_linshi <- p_linshi[which(is.na(p_linshi) == F)]
    if (length(p_linshi) > 1){
      P_PT_LNM_full[i,j] <- sumlog(p = p_linshi)$p
    }
    else if (length(p_linshi) == 1){
      P_PT_LNM_full[i,j] <- p_linshi
    }
  }
}
for (i in 1:28){
  for (j in 1:28){
    p_linshi <- c()
    for (q in 1:length(p_LNMT)){
      linshi <- p_LNMT[[q]]
      linshi[which(linshi == "NaN")] <- NA
      linshi1 <- matrix(NA, nrow=length(cell), ncol=length(cell),
                        dimnames=list(cell, cell))
      linshi1[rownames(linshi), colnames(linshi)] <- as.matrix(linshi)
      p_linshi <- c(p_linshi,
                    linshi1[rownames(P_LNMT_full)[i],colnames(P_LNMT_full)[j]])
    }
    p_linshi <- p_linshi[which(is.na(p_linshi) == F)]
    if (length(p_linshi) > 1){
      P_LNMT_full[i,j] <- sumlog(p = p_linshi)$p
    }
    else if (length(p_linshi) == 1){
      P_LNMT_full[i,j] <- p_linshi
    }
  }
}
P_PT_full1 <- P_PT_full
P_PT_LNM_full1 <- P_PT_LNM_full
P_LNMT_full1 <- P_LNMT_full
P_PT_full1[score_PT_full > 0] <- 1
P_PT_full1[score_PT_full < 0] <- -1
P_PT_LNM_full1[score_PT_LNM_full > 0] <- 1
P_PT_LNM_full1[score_PT_LNM_full < 0] <- -1
P_LNMT_full1[score_LNMT_full > 0] <- 1
P_LNMT_full1[score_LNMT_full < 0] <- -1
name <- c("Malignant C1","Malignant C2","Malignant C3","Malignant C4",
          "Malignant C5","Malignant C6","Malignant C7","Malignant C8",
          "Malignant C9","Malignant C10","Malignant C11","Malignant C12",
          "Malignant C13","Fibroblast","Endothelial",
          "Macrophages","Effector CD4+ T cell",
          "Effector CD8+ T cell","Inflammatory Treg",
          "Exhausted Treg","Exhausted CD8+ T cell",
          "DNT","Memory B cell",
          "Naive B cell","Proliferative B cell",
          "Transitional B cell","Plasma cell","Undefined")
score_PT_full <- score_PT_full[name,name]
score_PT_LNM_full <- score_PT_LNM_full[name,name]
score_LNMT_full <- score_LNMT_full[name,name]
P_PT_full <- P_PT_full[name,name]
P_PT_full1 <- P_PT_full1[name,name]
P_PT_LNM_full <- P_PT_LNM_full[name,name]
P_PT_LNM_full1 <- P_PT_LNM_full1[name,name]
P_LNMT_full <- P_LNMT_full[name,name]
P_LNMT_full1 <- P_LNMT_full1[name,name]
# PT heatmap
pheatmap(score_PT_full,fontsize=6,cellwidth = 12,cellheight = 8,
         color  = colorRampPalette(c("#92A8D1","white","#FF6347"))(100),
         border_color = "grey60",
         breaks = seq(-5, 5, length.out = 101),
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
col_fun <- colorRamp2(c(-1,0,1), c("#92A8D1", "white", "#FF6347"))
P_PT_full <- as.matrix(P_PT_full)
P_PT_full1 <- as.matrix(P_PT_full1)
P_PT_full1[P_PT_full > 0.05] <- 0
Heatmap(P_PT_full1,
        name = "P_PT_full1",
        col = col_fun,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        width = unit(12, "mm") * ncol(P_PT_full),
        height = unit(7, "mm") * nrow(P_PT_full),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.3f", P_PT_full[i, j]),
                    x, y, gp = gpar(fontsize = 10))
        })
# PT-LNM heatmap
pheatmap(score_PT_LNM_full,fontsize=6,cellwidth = 12,cellheight = 8,
         color  = colorRampPalette(c("#92A8D1","white","#FF6347"))(100),
         border_color = "grey60",
         breaks = seq(-5, 5, length.out = 101),
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
col_fun <- colorRamp2(c(-1,0,1), c("#92A8D1", "white", "#FF6347"))
P_PT_LNM_full <- as.matrix(P_PT_LNM_full)
P_PT_LNM_full1 <- as.matrix(P_PT_LNM_full1)
P_PT_LNM_full1[P_PT_LNM_full > 0.05] <- 0
Heatmap(P_PT_LNM_full1,
        name = "P_PT_LNM_full1",
        col = col_fun,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        width = unit(12, "mm") * ncol(P_PT_LNM_full),
        height = unit(7, "mm") * nrow(P_PT_LNM_full),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.3f", P_PT_LNM_full[i, j]),
                    x, y, gp = gpar(fontsize = 10))
        })
# LNMT heatmap
pheatmap(score_LNMT_full,fontsize=6,cellwidth = 12,cellheight = 8,
         color  = colorRampPalette(c("#92A8D1","white","#FF6347"))(100),
         border_color = "grey60",
         breaks = seq(-5, 5, length.out = 101),
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
col_fun <- colorRamp2(c(-1,0,1), c("#92A8D1", "white", "#FF6347"))
P_LNMT_full <- as.matrix(P_LNMT_full)
P_LNMT_full1 <- as.matrix(P_LNMT_full1)
P_LNMT_full1[P_LNMT_full > 0.05] <- 0
Heatmap(P_LNMT_full1,
        name = "P_LNMT_full1",
        col = col_fun,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        width = unit(12, "mm") * ncol(P_LNMT_full),
        height = unit(7, "mm") * nrow(P_LNMT_full),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.3f", P_LNMT_full[i, j]),
                    x, y, gp = gpar(fontsize = 10))
        })


#  Figure 3C
library(Seurat)
library(Hmisc)
library(pheatmap)
PT <- readRDS("PT_harmony8_afterQC_Doubcell_result.RDS")
PT_LNM <- readRDS("PT_LNM_harmony8_afterQC_Doubcell_result.RDS")
LNMT <- readRDS("LNMT_harmony8_afterQC_Doubcell_result.RDS")
PT_average <- table(PT$TMA,PT$cluster)
PT_LNM_average <- table(PT_LNM$TMA,PT_LNM$cluster)
LNMT_average <- table(LNMT$TMA,LNMT$cluster)
name <- c("Malignant C1","Malignant C2","Malignant C3","Malignant C4",
          "Malignant C5","Malignant C6","Malignant C7","Malignant C8",
          "Malignant C9","Malignant C10","Malignant C11","Malignant C12",
          "Malignant C13","Fibroblast","Endothelial",
          "Macrophages","Effector CD4+ T cell",
          "Effector CD8+ T cell","Inflammatory Treg",
          "Exhausted Treg","Exhausted CD8+ T cell",
          "DNT","Memory B cell",
          "Naive B cell","Proliferative B cell",
          "Transitional B cell","Plasma cell","Undefined")
#  calculation
PT_average <- as.matrix(PT_average)
PT_LNM_average <- as.matrix(PT_LNM_average)
LNMT_average <- as.matrix(LNMT_average)
PT_average <- PT_average[,name]
PT_LNM_average <- PT_LNM_average[,name]
LNMT_average <- LNMT_average[,name]
PT_cor <- rcorr(as.matrix(PT_average),type = "spearman")
PT_cor_r <- PT_cor$r
PT_cor_P <- PT_cor$P
PT_LNM_cor <- rcorr(as.matrix(PT_LNM_average),type = "spearman")
PT_LNM_cor_r <- PT_LNM_cor$r
PT_LNM_cor_P <- PT_LNM_cor$P
LNMT_cor <- rcorr(as.matrix(LNMT_average),type = "spearman")
LNMT_cor_r <- LNMT_cor$r
LNMT_cor_P <- LNMT_cor$P
# PT heatmap
pheatmap(PT_cor_r,fontsize=6,cellwidth = 12,cellheight = 8,
         color  = colorRampPalette(c("#92A8D1","white","#FF6347"))(100),
         border_color = "grey60",
         breaks = seq(-1, 1, length.out = 101),
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
# PT-LNM heatmap
pheatmap(PT_LNM_cor_r,fontsize=6,cellwidth = 12,cellheight = 8,
         color  = colorRampPalette(c("#92A8D1","white","#FF6347"))(100),
         border_color = "grey60",
         breaks = seq(-1, 1, length.out = 101),
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
# LNMT heatmap
pheatmap(LNMT_cor_r,fontsize=6,cellwidth = 12,cellheight = 8,
         color  = colorRampPalette(c("#92A8D1","white","#FF6347"))(100),
         border_color = "grey60",
         breaks = seq(-1, 1, length.out = 101),
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)


# Figure 4D
library(Seurat)
library(tidyverse)
PT <- readRDS("PT_harmony8_afterQC_Doubcell_result.RDS")
PT_LNM <- readRDS("PT_LNM_harmony8_afterQC_Doubcell_result.RDS")
LNMT <- readRDS("LNMT_harmony8_afterQC_Doubcell_result.RDS")
# subcluster
PT_Macrophages <-  subset(PT,subset = cluster == "Macrophages")
PT_EffCD4 <-  subset(PT,subset = cluster == "Effector CD4+ T cell")
PT_EffCD8 <-  subset(PT,subset = cluster == "Effector CD8+ T cell")
PT_MemB <-  subset(PT,subset = cluster == "Memory B cell")
PT_Plasma <-  subset(PT,subset = cluster == "Plasma cell")
PT_LNM_Macrophages <-  subset(PT_LNM,subset = cluster == "Macrophages")
PT_LNM_EffCD4 <-  subset(PT_LNM,subset = cluster == "Effector CD4+ T cell")
PT_LNM_EffCD8 <-  subset(PT_LNM,subset = cluster == "Effector CD8+ T cell")
PT_LNM_MemB <-  subset(PT_LNM,subset = cluster == "Memory B cell")
PT_LNM_Plasma <-  subset(PT_LNM,subset = cluster == "Plasma cell")
LNMT_Macrophages <-  subset(LNMT,subset = cluster == "Macrophages")
LNMT_EffCD4 <-  subset(LNMT,subset = cluster == "Effector CD4+ T cell")
LNMT_EffCD8 <-  subset(LNMT,subset = cluster == "Effector CD8+ T cell")
LNMT_MemB <-  subset(LNMT,subset = cluster == "Memory B cell")
LNMT_Plasma <-  subset(LNMT,subset = cluster == "Plasma cell")
# Average expression
PT_Macrophages <- AverageExpression(PT_Macrophages,assays = "RNA", 
                                    return.seurat = TRUE,
                                    group.by = "TMA",slot = "data")
PT_Macrophages <- ScaleData(PT_Macrophages)
PT_Macrophages <- GetAssayData(PT_Macrophages, slot = "counts")
PT_EffCD4 <- AverageExpression(PT_EffCD4,assays = "RNA", 
                               return.seurat = TRUE,
                               group.by = "TMA",slot = "data")
PT_EffCD4 <- ScaleData(PT_EffCD4)
PT_EffCD4 <- GetAssayData(PT_EffCD4, slot = "counts")
PT_EffCD8 <- AverageExpression(PT_EffCD8,assays = "RNA", 
                               return.seurat = TRUE,
                               group.by = "TMA",slot = "data")
PT_EffCD8 <- ScaleData(PT_EffCD8)
PT_EffCD8 <- GetAssayData(PT_EffCD8, slot = "counts")
PT_MemB <- AverageExpression(PT_MemB,assays = "RNA", 
                             return.seurat = TRUE,
                             group.by = "TMA",slot = "data")
PT_MemB <- ScaleData(PT_MemB)
PT_MemB <- GetAssayData(PT_MemB, slot = "counts")
PT_Plasma <- AverageExpression(PT_Plasma,assays = "RNA", 
                               return.seurat = TRUE,
                               group.by = "TMA",slot = "data")
PT_Plasma <- ScaleData(PT_Plasma)
PT_Plasma <- GetAssayData(PT_Plasma, slot = "counts")

PT_LNM_Macrophages <- AverageExpression(PT_LNM_Macrophages,assays = "RNA", 
                                        return.seurat = TRUE,
                                        group.by = "TMA",slot = "data")
PT_LNM_Macrophages <- ScaleData(PT_LNM_Macrophages)
PT_LNM_Macrophages <- GetAssayData(PT_LNM_Macrophages, slot = "counts")
PT_LNM_EffCD4 <- AverageExpression(PT_LNM_EffCD4,assays = "RNA", 
                                   return.seurat = TRUE,
                                   group.by = "TMA",slot = "data")
PT_LNM_EffCD4 <- ScaleData(PT_LNM_EffCD4)
PT_LNM_EffCD4 <- GetAssayData(PT_LNM_EffCD4, slot = "counts")
PT_LNM_EffCD8 <- AverageExpression(PT_LNM_EffCD8,assays = "RNA", 
                                   return.seurat = TRUE,
                                   group.by = "TMA",slot = "data")
PT_LNM_EffCD8 <- ScaleData(PT_LNM_EffCD8)
PT_LNM_EffCD8 <- GetAssayData(PT_LNM_EffCD8, slot = "counts")
PT_LNM_MemB <- AverageExpression(PT_LNM_MemB,assays = "RNA", 
                                 return.seurat = TRUE,
                                 group.by = "TMA",slot = "data")
PT_LNM_MemB <- ScaleData(PT_LNM_MemB)
PT_LNM_MemB <- GetAssayData(PT_LNM_MemB, slot = "counts")
PT_LNM_Plasma <- AverageExpression(PT_LNM_Plasma,assays = "RNA", 
                                   return.seurat = TRUE,
                                   group.by = "TMA",slot = "data")
PT_LNM_Plasma <- ScaleData(PT_LNM_Plasma)
PT_LNM_Plasma <- GetAssayData(PT_LNM_Plasma, slot = "counts")

LNMT_Macrophages <- AverageExpression(LNMT_Macrophages,assays = "RNA", 
                                      return.seurat = TRUE,
                                      group.by = "TMA",slot = "data")
LNMT_Macrophages <- ScaleData(LNMT_Macrophages)
LNMT_Macrophages <- GetAssayData(LNMT_Macrophages, slot = "counts")
LNMT_EffCD4 <- AverageExpression(LNMT_EffCD4,assays = "RNA", 
                                 return.seurat = TRUE,
                                 group.by = "TMA",slot = "data")
LNMT_EffCD4 <- ScaleData(LNMT_EffCD4)
LNMT_EffCD4 <- GetAssayData(LNMT_EffCD4, slot = "counts")
LNMT_EffCD8 <- AverageExpression(LNMT_EffCD8,assays = "RNA", 
                                 return.seurat = TRUE,
                                 group.by = "TMA",slot = "data")
LNMT_EffCD8 <- ScaleData(LNMT_EffCD8)
LNMT_EffCD8 <- GetAssayData(LNMT_EffCD8, slot = "counts")
LNMT_MemB <- AverageExpression(LNMT_MemB,assays = "RNA", 
                               return.seurat = TRUE,
                               group.by = "TMA",slot = "data")
LNMT_MemB <- ScaleData(LNMT_MemB)
LNMT_MemB <- GetAssayData(LNMT_MemB, slot = "counts")
LNMT_Plasma <- AverageExpression(LNMT_Plasma,assays = "RNA", 
                                 return.seurat = TRUE,
                                 group.by = "TMA",slot = "data")
LNMT_Plasma <- ScaleData(LNMT_Plasma)
LNMT_Plasma <- GetAssayData(LNMT_Plasma, slot = "counts")
PT_tma_inter_avoid <- readRDS("PT_TMA.RDS")
PT_LNM_tma_inter_avoid <- readRDS("PT_LNM_TMA.RDS")
LNMT_tma_inter_avoid <- readRDS("LNMT_TMA.RDS")
PT_tma_inter_avoid$TMA <- paste("g",PT_tma_inter_avoid$TMA,sep = "")
PT_LNM_tma_inter_avoid$TMA <- paste("g",PT_LNM_tma_inter_avoid$TMA,sep = "")
LNMT_tma_inter_avoid$TMA <- paste("g",LNMT_tma_inter_avoid$TMA,sep = "")
# PT
PT_Endothelial_Macrophages <- PT_tma_inter_avoid[which(PT_tma_inter_avoid$inter_cell == "Endothelial_Macrophages"),]
PT_Endothelial_effCD4 <- PT_tma_inter_avoid[which(PT_tma_inter_avoid$inter_cell == "Endothelial_effCD4"),]
PT_Endothelial_effCD8 <- PT_tma_inter_avoid[which(PT_tma_inter_avoid$inter_cell == "Endothelial_effCD8"),]
PT_Endothelial_MemB <- PT_tma_inter_avoid[which(PT_tma_inter_avoid$inter_cell == "Endothelial_MemB"),]
PT_Endothelial_Plasma <- PT_tma_inter_avoid[which(PT_tma_inter_avoid$inter_cell == "Endothelial_Plasma"),]
PT_Endothelial_Macrophages$CD68 <- PT_Macrophages["CD68",PT_Endothelial_Macrophages$TMA]
PT_Endothelial_Macrophages$CD163 <- PT_Macrophages["CD163",PT_Endothelial_Macrophages$TMA]
PT_Endothelial_effCD4$GZMA <- PT_EffCD4["GZMA",PT_Endothelial_effCD4$TMA]
PT_Endothelial_effCD4$GZMB <- PT_EffCD4["GZMB",PT_Endothelial_effCD4$TMA]
PT_Endothelial_effCD4$TCF7 <- PT_EffCD4["TCF7",PT_Endothelial_effCD4$TMA]
PT_Endothelial_effCD4$PDCD1 <- PT_EffCD4["PDCD1",PT_Endothelial_effCD4$TMA]
PT_Endothelial_effCD4$CTLA4 <- PT_EffCD4["CTLA4",PT_Endothelial_effCD4$TMA]
PT_Endothelial_effCD8$GZMA <- PT_EffCD8["GZMA",PT_Endothelial_effCD8$TMA]
PT_Endothelial_effCD8$GZMB <- PT_EffCD8["GZMB",PT_Endothelial_effCD8$TMA]
PT_Endothelial_effCD8$TCF7 <- PT_EffCD8["TCF7",PT_Endothelial_effCD8$TMA]
PT_Endothelial_effCD8$PDCD1 <- PT_EffCD8["PDCD1",PT_Endothelial_effCD8$TMA]
PT_Endothelial_effCD8$CTLA4 <- PT_EffCD8["CTLA4",PT_Endothelial_effCD8$TMA]
PT_Endothelial_MemB$CD19 <- PT_MemB["CD19",PT_Endothelial_MemB$TMA]
PT_Endothelial_MemB$CD27 <- PT_MemB["CD27",PT_Endothelial_MemB$TMA]
PT_Endothelial_MemB$CD80 <- PT_MemB["CD80",PT_Endothelial_MemB$TMA]
PT_Endothelial_Plasma$CD38 <- PT_Plasma["CD38",PT_Endothelial_Plasma$TMA]
# PT-LNM
PT_LNM_Endothelial_Macrophages <- PT_LNM_tma_inter_avoid[which(PT_LNM_tma_inter_avoid$inter_cell == "Endothelial_Macrophages"),]
PT_LNM_Endothelial_effCD4 <- PT_LNM_tma_inter_avoid[which(PT_LNM_tma_inter_avoid$inter_cell == "Endothelial_effCD4"),]
PT_LNM_Endothelial_effCD8 <- PT_LNM_tma_inter_avoid[which(PT_LNM_tma_inter_avoid$inter_cell == "Endothelial_effCD8"),]
PT_LNM_Endothelial_MemB <- PT_LNM_tma_inter_avoid[which(PT_LNM_tma_inter_avoid$inter_cell == "Endothelial_MemB"),]
PT_LNM_Endothelial_Plasma <- PT_LNM_tma_inter_avoid[which(PT_LNM_tma_inter_avoid$inter_cell == "Endothelial_Plasma"),]
PT_LNM_Endothelial_Macrophages$CD68 <- PT_LNM_Macrophages["CD68",PT_LNM_Endothelial_Macrophages$TMA]
PT_LNM_Endothelial_Macrophages$CD163 <- PT_LNM_Macrophages["CD163",PT_LNM_Endothelial_Macrophages$TMA]
PT_LNM_Endothelial_effCD4$GZMA <- PT_LNM_EffCD4["GZMA",PT_LNM_Endothelial_effCD4$TMA]
PT_LNM_Endothelial_effCD4$GZMB <- PT_LNM_EffCD4["GZMB",PT_LNM_Endothelial_effCD4$TMA]
PT_LNM_Endothelial_effCD4$TCF7 <- PT_LNM_EffCD4["TCF7",PT_LNM_Endothelial_effCD4$TMA]
PT_LNM_Endothelial_effCD4$PDCD1 <- PT_LNM_EffCD4["PDCD1",PT_LNM_Endothelial_effCD4$TMA]
PT_LNM_Endothelial_effCD4$CTLA4 <- PT_LNM_EffCD4["CTLA4",PT_LNM_Endothelial_effCD4$TMA]
PT_LNM_Endothelial_effCD8$GZMA <- PT_LNM_EffCD8["GZMA",PT_LNM_Endothelial_effCD8$TMA]
PT_LNM_Endothelial_effCD8$GZMB <- PT_LNM_EffCD8["GZMB",PT_LNM_Endothelial_effCD8$TMA]
PT_LNM_Endothelial_effCD8$TCF7 <- PT_LNM_EffCD8["TCF7",PT_LNM_Endothelial_effCD8$TMA]
PT_LNM_Endothelial_effCD8$PDCD1 <- PT_LNM_EffCD8["PDCD1",PT_LNM_Endothelial_effCD8$TMA]
PT_LNM_Endothelial_effCD8$CTLA4 <- PT_LNM_EffCD8["CTLA4",PT_LNM_Endothelial_effCD8$TMA]
PT_LNM_Endothelial_MemB$CD19 <- PT_LNM_MemB["CD19",PT_LNM_Endothelial_MemB$TMA]
PT_LNM_Endothelial_MemB$CD27 <- PT_LNM_MemB["CD27",PT_LNM_Endothelial_MemB$TMA]
PT_LNM_Endothelial_MemB$CD80 <- PT_LNM_MemB["CD80",PT_LNM_Endothelial_MemB$TMA]
PT_LNM_Endothelial_Plasma$CD38 <- PT_LNM_Plasma["CD38",PT_LNM_Endothelial_Plasma$TMA]
# LNMT
LNMT_Endothelial_Macrophages <- LNMT_tma_inter_avoid[which(LNMT_tma_inter_avoid$inter_cell == "Endothelial_Macrophages"),]
LNMT_Endothelial_effCD4 <- LNMT_tma_inter_avoid[which(LNMT_tma_inter_avoid$inter_cell == "Endothelial_effCD4"),]
LNMT_Endothelial_effCD8 <- LNMT_tma_inter_avoid[which(LNMT_tma_inter_avoid$inter_cell == "Endothelial_effCD8"),]
LNMT_Endothelial_MemB <- LNMT_tma_inter_avoid[which(LNMT_tma_inter_avoid$inter_cell == "Endothelial_MemB"),]
LNMT_Endothelial_Plasma <- LNMT_tma_inter_avoid[which(LNMT_tma_inter_avoid$inter_cell == "Endothelial_Plasma"),]
LNMT_Endothelial_Macrophages$CD68 <- LNMT_Macrophages["CD68",LNMT_Endothelial_Macrophages$TMA]
LNMT_Endothelial_Macrophages$CD163 <- LNMT_Macrophages["CD163",LNMT_Endothelial_Macrophages$TMA]
LNMT_Endothelial_effCD4$GZMA <- LNMT_EffCD4["GZMA",LNMT_Endothelial_effCD4$TMA]
LNMT_Endothelial_effCD4$GZMB <- LNMT_EffCD4["GZMB",LNMT_Endothelial_effCD4$TMA]
LNMT_Endothelial_effCD4$TCF7 <- LNMT_EffCD4["TCF7",LNMT_Endothelial_effCD4$TMA]
LNMT_Endothelial_effCD4$PDCD1 <- LNMT_EffCD4["PDCD1",LNMT_Endothelial_effCD4$TMA]
LNMT_Endothelial_effCD4$CTLA4 <- LNMT_EffCD4["CTLA4",LNMT_Endothelial_effCD4$TMA]
LNMT_Endothelial_effCD8$GZMA <- LNMT_EffCD8["GZMA",LNMT_Endothelial_effCD8$TMA]
LNMT_Endothelial_effCD8$GZMB <- LNMT_EffCD8["GZMB",LNMT_Endothelial_effCD8$TMA]
LNMT_Endothelial_effCD8$TCF7 <- LNMT_EffCD8["TCF7",LNMT_Endothelial_effCD8$TMA]
LNMT_Endothelial_effCD8$PDCD1 <- LNMT_EffCD8["PDCD1",LNMT_Endothelial_effCD8$TMA]
LNMT_Endothelial_effCD8$CTLA4 <- LNMT_EffCD8["CTLA4",LNMT_Endothelial_effCD8$TMA]
LNMT_Endothelial_MemB$CD19 <- LNMT_MemB["CD19",LNMT_Endothelial_MemB$TMA]
LNMT_Endothelial_MemB$CD27 <- LNMT_MemB["CD27",LNMT_Endothelial_MemB$TMA]
LNMT_Endothelial_MemB$CD80 <- LNMT_MemB["CD80",LNMT_Endothelial_MemB$TMA]
LNMT_Endothelial_Plasma$CD38 <- LNMT_Plasma["CD38",LNMT_Endothelial_Plasma$TMA]
##  function defined
library(ggplot2)
library(gghalves)
library(ggpubr)
com <- list( c("PT_Interact", "PT_Avoid"), 
             c("PT_LNM_Interact", "PT_LNM_Avoid"), 
             c("LNMT_Interact", "LNMT_Avoid") )
options(warn = -1)
ALL_create_half_boxplot <- function(data, x, y, fill_colors) {
  ggplot(data, aes_string(x = x, y = paste0("log2(", y, "+ 1)"), fill = x)) + 
    geom_half_boxplot(width = 0.4, position = position_dodge(0.7), alpha = 1, outlier.shape = NA) +
    geom_half_point(aes_string(color = x), shape = 16, alpha = 1, size = 3,
                    transformation = position_jitter(width = 0.1,height = 0)) +
    scale_fill_manual(values = c("#E74C3C","#3498DB","#E74C3C","#3498DB","#E74C3C","#3498DB")) +
    scale_color_manual(values = c("#E74C3C","#3498DB","#E74C3C","#3498DB","#E74C3C","#3498DB")) +
    stat_compare_means(comparisons = com,method = "wilcox.test") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.text.y = element_text(angle = 0, hjust = 1),
          legend.position = "right") +
    labs(x = '', y = y)
}
###  Endothelial_effCD4 
Endothelial_effCD4 <- as.data.frame(rbind(PT_Endothelial_effCD4,
                                          PT_LNM_Endothelial_effCD4,
                                          LNMT_Endothelial_effCD4))
Endothelial_effCD4$mer_type <- paste(Endothelial_effCD4$meta,
                                     Endothelial_effCD4$inter_type,
                                     sep = "_")
Endothelial_effCD4$mer_type <- factor(Endothelial_effCD4$mer_type,
                                      levels = c("PT_Interact", "PT_Avoid",
                                                 "PT_LNM_Interact", "PT_LNM_Avoid",
                                                 "LNMT_Interact", "LNMT_Avoid"))
ALL_create_half_boxplot(Endothelial_effCD4,"mer_type","CTLA4","inter_type")
ALL_create_half_boxplot(Endothelial_effCD4,"mer_type","GZMA","inter_type")
ALL_create_half_boxplot(Endothelial_effCD4,"mer_type","GZMB","inter_type")
###  Endothelial_effCD8 
Endothelial_effCD8 <- as.data.frame(rbind(PT_Endothelial_effCD8,
                                          PT_LNM_Endothelial_effCD8,
                                          LNMT_Endothelial_effCD8))
Endothelial_effCD8$mer_type <- paste(Endothelial_effCD8$meta,
                                     Endothelial_effCD8$inter_type,
                                     sep = "_")
Endothelial_effCD8$mer_type <- factor(Endothelial_effCD8$mer_type,
                                      levels = c("PT_Interact", "PT_Avoid",
                                                 "PT_LNM_Interact", "PT_LNM_Avoid",
                                                 "LNMT_Interact", "LNMT_Avoid"))
ALL_create_half_boxplot(Endothelial_effCD8,"mer_type","PDCD1","inter_type")
ALL_create_half_boxplot(Endothelial_effCD8,"mer_type","CTLA4","inter_type")
ALL_create_half_boxplot(Endothelial_effCD8,"mer_type","GZMA","inter_type")
ALL_create_half_boxplot(Endothelial_effCD8,"mer_type","GZMB","inter_type")
###  Endothelial_Macrophages
Endothelial_Macrophages <- as.data.frame(rbind(PT_Endothelial_Macrophages,
                                               PT_LNM_Endothelial_Macrophages,
                                               LNMT_Endothelial_Macrophages))
Endothelial_Macrophages$mer_type <- paste(Endothelial_Macrophages$meta,
                                          Endothelial_Macrophages$inter_type,
                                          sep = "_")
Endothelial_Macrophages$mer_type <- factor(Endothelial_Macrophages$mer_type,
                                           levels = c("PT_Interact", "PT_Avoid",
                                                      "PT_LNM_Interact", "PT_LNM_Avoid",
                                                      "LNMT_Interact", "LNMT_Avoid"))
ALL_create_half_boxplot(Endothelial_Macrophages,"mer_type","CD68","inter_type")
ALL_create_half_boxplot(Endothelial_Macrophages,"mer_type","CD163","inter_type")
###  Endothelial_MemB
Endothelial_MemB <- as.data.frame(rbind(PT_Endothelial_MemB,
                                        PT_LNM_Endothelial_MemB,
                                        LNMT_Endothelial_MemB))
Endothelial_MemB$mer_type <- paste(Endothelial_MemB$meta,
                                   Endothelial_MemB$inter_type,
                                   sep = "_")
Endothelial_MemB$mer_type <- factor(Endothelial_MemB$mer_type,
                                    levels = c("PT_Interact", "PT_Avoid",
                                               "PT_LNM_Interact", "PT_LNM_Avoid",
                                               "LNMT_Interact", "LNMT_Avoid"))
ALL_create_half_boxplot(Endothelial_MemB,"mer_type","CD19","inter_type")
ALL_create_half_boxplot(Endothelial_MemB,"mer_type","CD27","inter_type")
ALL_create_half_boxplot(Endothelial_MemB,"mer_type","CD80","inter_type")
###  Endothelial_Plasma
Endothelial_Plasma <- as.data.frame(rbind(PT_Endothelial_Plasma,
                                          PT_LNM_Endothelial_Plasma,
                                          LNMT_Endothelial_Plasma))
Endothelial_Plasma$mer_type <- paste(Endothelial_Plasma$meta,
                                     Endothelial_Plasma$inter_type,
                                     sep = "_")
Endothelial_Plasma$mer_type <- factor(Endothelial_Plasma$mer_type,
                                      levels = c("PT_Interact", "PT_Avoid",
                                                 "PT_LNM_Interact", "PT_LNM_Avoid",
                                                 "LNMT_Interact", "LNMT_Avoid"))
ALL_create_half_boxplot(Endothelial_Plasma,"mer_type","CD38","inter_type")



















