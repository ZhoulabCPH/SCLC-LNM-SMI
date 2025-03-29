
# Figure 3A
library(Seurat)
library(tidyverse)
library(ggplot2)
library(harmony)
DimPlot(immuneT, reduction = "umap", 
        cols = c("#BDA7CB","#ABCADE","#3B76AD",
                 "#C0937E","#DA8F6F","#67A59B"),
        group.by = "clusters",raster=T,pt.size = 1.5) 
immuneT$type <- paste(immuneT$Patient_metastasis,
                      immuneT$Tissue_metastasis)
immuneT_PT <- subset(immuneT,subset = type == "Primary Primary")
immuneT_PT_LNM <- subset(immuneT,subset = type == "Metastasis Primary")
immuneT_LNMT <- subset(immuneT,subset = type == "Metastasis LNM")
DimPlot(immuneT_PT, reduction = "umap", 
        cols = c("#BDA7CB","#ABCADE","#3B76AD",
                 "#C0937E","#DA8F6F","#67A59B"),
        group.by = "clusters",raster=T,pt.size = 1.5) 
DimPlot(immuneT_PT_LNM, reduction = "umap", 
        cols = c("#BDA7CB","#ABCADE","#3B76AD",
                 "#C0937E","#DA8F6F","#67A59B"),
        group.by = "clusters",raster=T,pt.size = 1.5) 
DimPlot(immuneT_LNMT, reduction = "umap", 
        cols = c("#BDA7CB","#ABCADE","#3B76AD",
                 "#C0937E","#DA8F6F","#67A59B"),
        group.by = "clusters",raster=T,pt.size = 1.5) 
markers <- c("CD4","CD8A", "CD8B","FOXP3",
             "IFNGR1","CXCL10","CXCL12",
             "CCL21","CXCL14","HSPA1A/B",
             "GZMA","GZMB","GZMK","IFNG",
             "PDCD1","HAVCR2","TIGIT","TOX")
immuneT_data <- AverageExpression(immuneT,
                                  assays = "RNA", 
                                  return.seurat = TRUE,
                                  group.by = "clusters",
                                  slot = "data")
cluster.averages <- ScaleData(immuneT_data)
mat <- GetAssayData(cluster.averages, slot = "scale.data")
mat <- as.matrix(mat[markers, ])
mat1 <- GetAssayData(immuneT_data, slot = "data")
mat1 <- as.matrix(mat1[c("CD3D","CD3E"), ])
mat[mat > 1] <- 1
mat[mat < (0)] <- (0)
library(pheatmap)
pheatmap(t(mat),fontsize=6,gaps_col = c(4,10,14),
         cellwidth = 26,cellheight = 26,
         color  = colorRampPalette(c("white","#483D8B"))(100),
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
breaks <- seq(0, 4, length.out = 101)
pheatmap(t(mat1),fontsize=6,cellwidth = 26,cellheight = 26,
         color  = colorRampPalette(c("white","#8B0000"))(100),
         border_color = "grey60",breaks = breaks,
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)

# Figure 3B 
library(Seurat)
library(tidyverse)
library(ggplot2)
library(harmony)
DimPlot(immuneB, reduction = "umap", 
        cols = c("#58888C","#74A764",
                 "#F5B8A8","#A82D06"),
        group.by = "clusters",raster=T,pt.size = 2) 
immuneB$type <- paste(immuneB$Patient_metastasis,
                      immuneB$Tissue_metastasis)
immuneB_PT <- subset(immuneB,subset = type == "Primary Primary")
immuneB_PT_LNM <- subset(immuneB,subset = type == "Metastasis Primary")
immuneB_LNMT <- subset(immuneB,subset = type == "Metastasis LNM")
DimPlot(immuneB_PT, reduction = "umap", 
        cols = c("#58888C","#74A764",
                 "#F5B8A8","#A82D06"),
        group.by = "clusters",raster=T,pt.size = 3) 
DimPlot(immuneB_PT_LNM, reduction = "umap", 
        cols = c("#58888C","#74A764",
                 "#F5B8A8","#A82D06"),
        group.by = "clusters",raster=T,pt.size = 3) 
DimPlot(immuneB_LNMT, reduction = "umap", 
        cols = c("#58888C","#74A764",
                 "#F5B8A8","#A82D06"),
        group.by = "clusters",raster=T,pt.size = 3) 
markers <- c("MKI67","STMN1","PCNA",  # Proliferative B cell
             "SELL","CD83","CD22","TCL1A", # Naive B cell
             "CTSD","IGHA1","IGHG2") # Memory B cell
immuneB_data <- AverageExpression(immuneB,
                                    assays = "RNA", 
                                    return.seurat = TRUE,
                                    group.by = "clusters",
                                    slot = "data")
cluster.averages <- ScaleData(immuneB_data)
mat <- GetAssayData(cluster.averages, slot = "scale.data")
mat <- as.matrix(mat[markers, ])
mat[mat > 1] <- 1
mat[mat < (0)] <- (0)
library(pheatmap)
pheatmap(t(mat),fontsize=6,#gaps_col = c(3,6,9),
         cellwidth = 26,cellheight = 26,
         color  = colorRampPalette(c("white","#483D8B"))(100),
         border_color = "grey60",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)


# Figure 3C
library(Seurat)
library(tidyverse)
library(ggplot2)
Tcell$type <- paste(Tcell$Patient_metastasis,Tcell$Tissue_metastasis)
Bcell$type <- paste(Bcell$Patient_metastasis,Bcell$Tissue_metastasis)
# T cell
Tcell_PT <- subset(Tcell,subset = type == "Primary Primary")
Tcell_PT_LNM <- subset(Tcell,subset = type == "Metastasis Primary")
Tcell_LNMT <- subset(Tcell,subset = type == "Metastasis LNM")
# B cell
Bcell_PT <- subset(Bcell,subset = type == "Primary Primary")
Bcell_PT_LNM <- subset(Bcell,subset = type == "Metastasis Primary")
Bcell_LNMT <- subset(Bcell,subset = type == "Metastasis LNM")
#  PT, PT-LNM and LNMT  T cell
pt_rate <- table(Tcell_PT$clusters) / sum(table(Tcell_PT$clusters))
ptlnm_rate <- table(Tcell_PT_LNM$clusters) / sum(table(Tcell_PT_LNM$clusters))
lnmt_rate <- table(Tcell_LNMT$clusters) / sum(table(Tcell_LNMT$clusters))
total <- data.frame(rate = c(pt_rate,ptlnm_rate,lnmt_rate),
                    type = c(rep("PT",length(pt_rate)),
                             rep("PT-LNM",length(ptlnm_rate)),
                             rep("LNMT",length(lnmt_rate))),
                    cell = c(names(pt_rate),names(ptlnm_rate),
                             names(lnmt_rate)))
total$cell <- factor(total$cell,levels = names(pt_rate))
total$type <- factor(total$type,levels = c("PT","PT-LNM","LNMT"))
ggplot(total, aes(x = type, y = rate, fill = cell)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#BDA7CB","#ABCADE","#3B76AD",
                               "#C0937E","#DA8F6F","#67A59B"))+
  theme_bw() + ylim(0,1)+
  theme(axis.text.x = element_text(hjust = 1, colour = "black", size = 10, angle = 90), 
        axis.text.y = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10, colour = "black"), 
        legend.text = element_text(colour = "black", size = 10), 
        legend.title = element_text(colour = "black", size = 10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() ) +
  ylab("Percentage (%)")
#  PT, PT-LNM and LNMT B cell 
pt_rate <- table(Bcell_PT$clusters) / sum(table(Bcell_PT$clusters))
ptlnm_rate <- table(Bcell_PT_LNM$clusters) / sum(table(Bcell_PT_LNM$clusters))
lnmt_rate <- table(Bcell_LNMT$clusters) / sum(table(Bcell_LNMT$clusters))
total <- data.frame(rate = c(pt_rate,ptlnm_rate,lnmt_rate),
                    type = c(rep("PT",length(pt_rate)),
                             rep("PT-LNM",length(ptlnm_rate)),
                             rep("LNMT",length(lnmt_rate))),
                    cell = c(names(pt_rate),names(ptlnm_rate),
                             names(lnmt_rate)))
total$cell <- factor(total$cell,levels = names(pt_rate))
total$type <- factor(total$type,levels = c("PT","PT-LNM","LNMT"))
ggplot(total, aes(x = type, y = rate, fill = cell)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#58888C","#74A764",
                               "#F5B8A8","#A82D06"))+
  theme_bw() + ylim(0,1)+
  theme(axis.text.x = element_text(hjust = 1, colour = "black", size = 10, angle = 90), 
        axis.text.y = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10, colour = "black"), 
        legend.text = element_text(colour = "black", size = 10), 
        legend.title = element_text(colour = "black", size = 10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() ) +
  ylab("Percentage (%)")


#  Figure 3D
library(Seurat)
library(tidyverse)
library(ggplot2)
Tcell <- readRDS("Tcell2d5_harmony8_afterQC_Doubcell_result.RDS")
Bcell <- readRDS("Bcell2d5_harmony8_afterQC_Doubcell_result.RDS")
Tcell$type <- paste(Tcell$Patient_metastasis,
                    Tcell$Tissue_metastasis)
Bcell$type <- paste(Bcell$Patient_metastasis,
                    Bcell$Tissue_metastasis)
# T cell
Tcell_PT <- subset(Tcell,subset = type == "Primary Primary")
Tcell_PT_LNM <- subset(Tcell,subset = type == "Metastasis Primary")
Tcell_LNMT <- subset(Tcell,subset = type == "Metastasis LNM")
# B cell
Bcell_PT <- subset(Bcell,subset = type == "Primary Primary")
Bcell_PT_LNM <- subset(Bcell,subset = type == "Metastasis Primary")
Bcell_LNMT <- subset(Bcell,subset = type == "Metastasis LNM")
#  PT, PT-LNM 和 LNMT T cell 
pt_rateT <- table(Tcell_PT$clusters) / sum(table(Tcell_PT$clusters))
ptlnm_rateT <- table(Tcell_PT_LNM$clusters) / sum(table(Tcell_PT_LNM$clusters))
lnmt_rateT <- table(Tcell_LNMT$clusters) / sum(table(Tcell_LNMT$clusters))
#  PT, PT-LNM 和 LNMT B cell  
pt_rateB <- table(Bcell_PT$clusters) / sum(table(Bcell_PT$clusters))
ptlnm_rateB <- table(Bcell_PT_LNM$clusters) / sum(table(Bcell_PT_LNM$clusters))
lnmt_rateB <- table(Bcell_LNMT$clusters) / sum(table(Bcell_LNMT$clusters))
# FC 
PTLNM_vs_PT <- data.frame(log2fc = log2(c(((ptlnm_rateT) / (pt_rateT)),
                                          ((ptlnm_rateB) / (pt_rateB)))),
                          cell = c(names(ptlnm_rateT),
                                   names(ptlnm_rateB)))
LNMT_vs_PT <- data.frame(log2fc = log2(c(((lnmt_rateT) / (pt_rateT)),
                                         ((lnmt_rateB) / (pt_rateB)))),
                         cell = c(names(lnmt_rateT),
                                  names(lnmt_rateB)))
LNMT_vs_PTLNM <- data.frame(log2fc = log2(c(((lnmt_rateT) / (ptlnm_rateT)),
                                            ((lnmt_rateB) / (ptlnm_rateB)))),
                            cell = c(names(lnmt_rateT),
                                     names(lnmt_rateB)))
PTLNM_vs_PT$cell <- factor(PTLNM_vs_PT$cell,
                           levels = c("Effector CD4+ T cell","Effector CD8+ T cell",
                                      "Inflammatory Treg","Exhausted Treg",
                                      "Exhausted CD8+ T cell","DNT",
                                      "Memory B cell","Naive B cell",
                                      "Proliferative B cell","Transitional B cell"))
LNMT_vs_PTLNM$cell <- factor(LNMT_vs_PTLNM$cell,
                             levels = c("Effector CD4+ T cell","Effector CD8+ T cell",
                                        "Inflammatory Treg","Exhausted Treg",
                                        "Exhausted CD8+ T cell","DNT",
                                        "Memory B cell","Naive B cell",
                                        "Proliferative B cell","Transitional B cell"))
LNMT_vs_PT$cell <- factor(LNMT_vs_PT$cell,
                          levels = c("Effector CD4+ T cell","Effector CD8+ T cell",
                                     "Inflammatory Treg","Exhausted Treg",
                                     "Exhausted CD8+ T cell","DNT",
                                     "Memory B cell","Naive B cell",
                                     "Proliferative B cell","Transitional B cell"))
ggplot(data = PTLNM_vs_PT, aes(x = log2fc, y = cell, color = cell)) + 
  geom_point(size = 5, alpha = 0.7) + 
  xlim(-2, 5) + theme_bw() +
  scale_color_manual(values = c("DNT" = "#BDA7CB", 
                                "Effector CD4+ T cell" = "#ABCADE",
                                "Effector CD8+ T cell" = "#3B76AD",
                                "Exhausted CD8+ T cell" = "#C0937E",
                                "Exhausted Treg" = "#DA8F6F",
                                "Inflammatory Treg" = "#67A59B",
                                "Memory B cell" = "#58888C",
                                "Naive B cell" = "#74A764",
                                "Proliferative B cell" = "#F5B8A8",
                                "Transitional B cell" = "#A82D06"))
ggplot(data = LNMT_vs_PT, aes(x = log2fc, y = cell, color = cell)) + 
  geom_point(size = 5, alpha = 0.7) + 
  xlim(-2, 5) + theme_bw() +
  scale_color_manual(values = c("DNT" = "#BDA7CB", 
                                "Effector CD4+ T cell" = "#ABCADE",
                                "Effector CD8+ T cell" = "#3B76AD",
                                "Exhausted CD8+ T cell" = "#C0937E",
                                "Exhausted Treg" = "#DA8F6F",
                                "Inflammatory Treg" = "#67A59B",
                                "Memory B cell" = "#58888C",
                                "Naive B cell" = "#74A764",
                                "Proliferative B cell" = "#F5B8A8",
                                "Transitional B cell" = "#A82D06"))
ggplot(data = LNMT_vs_PTLNM, aes(x = log2fc, y = cell, color = cell)) + 
  geom_point(size = 5, alpha = 0.7) + 
  xlim(-2, 11) + theme_bw() +
  scale_color_manual(values = c("DNT" = "#BDA7CB", 
                                "Effector CD4+ T cell" = "#ABCADE",
                                "Effector CD8+ T cell" = "#3B76AD",
                                "Exhausted CD8+ T cell" = "#C0937E",
                                "Exhausted Treg" = "#DA8F6F",
                                "Inflammatory Treg" = "#67A59B",
                                "Memory B cell" = "#58888C",
                                "Naive B cell" = "#74A764",
                                "Proliferative B cell" = "#F5B8A8",
                                "Transitional B cell" = "#A82D06"))


# Figure 3E
library(Seurat)
library(tidyverse)
library(ggplot2)
PTLNMvsPT_deg <- data.frame(number = as.numeric(table(PTLNMvsPT_marker$cell)),
                            cell = names(table(PTLNMvsPT_marker$cell)))
LNMTvsPT_deg <- data.frame(number = as.numeric(table(LNMTvsPT_marker$cell)),
                           cell = names(table(LNMTvsPT_marker$cell)))
LNMTvsPTlnm_deg <- data.frame(number = as.numeric(table(LNMTvsPTlnm_marker$cell)),
                              cell = names(table(LNMTvsPTlnm_marker$cell)))

PTLNMvsPT_deg$cell <- factor(PTLNMvsPT_deg$cell,
                             levels = c("Effector CD4+ T cell","Effector CD8+ T cell",
                                        "Inflammatory Treg","Exhausted Treg",
                                        "Exhausted CD8+ T cell","DNT",
                                        "Memory B cell",
                                        "Proliferative B cell","Transitional B cell"))
LNMTvsPT_deg$cell <- factor(LNMTvsPT_deg$cell,
                            levels = c("Effector CD4+ T cell","Effector CD8+ T cell",
                                       "Inflammatory Treg","Exhausted Treg",
                                       "Exhausted CD8+ T cell","DNT",
                                       "Memory B cell",
                                       "Proliferative B cell","Transitional B cell"))
LNMTvsPTlnm_deg$cell <- factor(LNMTvsPTlnm_deg$cell,
                               levels = c("Effector CD4+ T cell","Effector CD8+ T cell",
                                          "Inflammatory Treg","Exhausted Treg",
                                          "Exhausted CD8+ T cell","DNT",
                                          "Memory B cell","Naive B cell",
                                          "Proliferative B cell","Transitional B cell"))
ggplot(data = PTLNMvsPT_deg, aes(x = number, y = cell, color = cell)) + 
  geom_point(size = 3, alpha = 0.7) + 
  xlim(0, max(PTLNMvsPT_deg$number)+2) + theme_bw() +
  scale_color_manual(values = c("DNT" = "#BDA7CB", 
                                "Effector CD4+ T cell" = "#ABCADE",
                                "Effector CD8+ T cell" = "#3B76AD",
                                "Exhausted CD8+ T cell" = "#C0937E",
                                "Exhausted Treg" = "#DA8F6F",
                                "Inflammatory Treg" = "#67A59B",
                                "Memory B cell" = "#58888C",
                                "Naive B cell" = "#74A764",
                                "Proliferative B cell" = "#F5B8A8",
                                "Transitional B cell" = "#A82D06"))
ggplot(data = LNMTvsPT_deg, aes(x = number, y = cell, color = cell)) + 
  geom_point(size = 3, alpha = 0.7) + 
  xlim(0, max(LNMTvsPT_deg$number)+2) + theme_bw() +
  scale_color_manual(values = c("DNT" = "#BDA7CB", 
                                "Effector CD4+ T cell" = "#ABCADE",
                                "Effector CD8+ T cell" = "#3B76AD",
                                "Exhausted CD8+ T cell" = "#C0937E",
                                "Exhausted Treg" = "#DA8F6F",
                                "Inflammatory Treg" = "#67A59B",
                                "Memory B cell" = "#58888C",
                                "Naive B cell" = "#74A764",
                                "Proliferative B cell" = "#F5B8A8",
                                "Transitional B cell" = "#A82D06"))
ggplot(data = LNMTvsPTlnm_deg, aes(x = number, y = cell, color = cell)) + 
  geom_point(size = 3, alpha = 0.7) + 
  xlim(0, max(LNMTvsPTlnm_deg$number)+2) + theme_bw() +
  scale_color_manual(values = c("DNT" = "#BDA7CB", 
                                "Effector CD4+ T cell" = "#ABCADE",
                                "Effector CD8+ T cell" = "#3B76AD",
                                "Exhausted CD8+ T cell" = "#C0937E",
                                "Exhausted Treg" = "#DA8F6F",
                                "Inflammatory Treg" = "#67A59B",
                                "Memory B cell" = "#58888C",
                                "Naive B cell" = "#74A764",
                                "Proliferative B cell" = "#F5B8A8",
                                "Transitional B cell" = "#A82D06"))

# Figure 3F
library(ggplot2)
library(ggstatsplot)
library(ggpubr)
library(ggExtra)
rownames(PTLNMvsPT_marker) <- paste(PTLNMvsPT_marker[,1],PTLNMvsPT_marker$cell)
rownames(LNMTvsPT_marker) <- paste(LNMTvsPT_marker[,1],LNMTvsPT_marker$cell)
rownames(LNMTvsPTlnm_marker) <- paste(LNMTvsPTlnm_marker[,1],LNMTvsPTlnm_marker$cell)
PTLNMvsPT_marker_up <- PTLNMvsPT_marker[which(PTLNMvsPT_marker$avg_log2FC > 0.5),]
PTLNMvsPT_marker_down <- PTLNMvsPT_marker[which(PTLNMvsPT_marker$avg_log2FC < 0.5),]
LNMTvsPT_marker_up <- LNMTvsPT_marker[which(LNMTvsPT_marker$avg_log2FC > 0.5),]
LNMTvsPT_marker_down <- LNMTvsPT_marker[which(LNMTvsPT_marker$avg_log2FC < 0.5),]
LNMTvsPTlnm_marker_up <- LNMTvsPTlnm_marker[which(LNMTvsPTlnm_marker$avg_log2FC > 0.5),]
LNMTvsPTlnm_marker_down <- LNMTvsPTlnm_marker[which(LNMTvsPTlnm_marker$avg_log2FC < 0.5),]
name1 <- intersect(PTLNMvsPT_marker_up[,1],LNMTvsPT_marker_up[,1])
name2 <- intersect(LNMTvsPTlnm_marker_up[,1],LNMTvsPT_marker_up[,1])
name3 <- intersect(LNMTvsPTlnm_marker_up[,1],PTLNMvsPT_marker_up[,1])
name11 <- intersect(PTLNMvsPT_marker_down[,1],LNMTvsPT_marker_down[,1])
name22 <- intersect(LNMTvsPTlnm_marker_down[,1],LNMTvsPT_marker_down[,1])
name33 <- intersect(LNMTvsPTlnm_marker_down[,1],PTLNMvsPT_marker_down[,1])
rownames(PTLNMvsPT_marker_up) <- PTLNMvsPT_marker_up[,1]
rownames(LNMTvsPT_marker_up) <- LNMTvsPT_marker_up[,1]
rownames(LNMTvsPTlnm_marker_up) <- LNMTvsPTlnm_marker_up[,1]
rownames(PTLNMvsPT_marker_down) <- PTLNMvsPT_marker_down[,1]
rownames(LNMTvsPT_marker_down) <- LNMTvsPT_marker_down[,1]
rownames(LNMTvsPTlnm_marker_down) <- LNMTvsPTlnm_marker_down[,1]
inter_mat1 <- data.frame(log2fc_PTLNMvsPT = c(PTLNMvsPT_marker_up[name1,3],
                                              PTLNMvsPT_marker_down[name11,3]),
                         log2fc_LNMTvsPT = c(LNMTvsPT_marker_up[name1,3],
                                             LNMTvsPT_marker_down[name11,3]))
inter_mat2 <- data.frame(log2fc_LNMTvsPTlnm = c(LNMTvsPTlnm_marker_up[name2,3],
                                                LNMTvsPTlnm_marker_down[name22,3]),
                         log2fc_LNMTvsPT = c(LNMTvsPT_marker_up[name2,3],
                                             LNMTvsPT_marker_down[name22,3]))
p1 <- ggplot(inter_mat1, aes(x=log2fc_PTLNMvsPT, y=log2fc_LNMTvsPT)) +       
  xlab("log2fc_PTLNMvsPT")+ylab("log2fc_LNMTvsPT")+   
  geom_point()+ geom_smooth(method="lm",colour='#764C29',fill='#E7E1D7') + 
  theme_bw()+      
  stat_cor(method = 'spearman', aes(x =log2fc_PTLNMvsPT, y =log2fc_LNMTvsPT))
ggMarginal(p1, type = "density", 
           xparams = list(fill = "#D3D3D3"),
           yparams = list(fill = "#D3D3D3"))
p2 <- ggplot(inter_mat2, aes(x=log2fc_LNMTvsPTlnm, y=log2fc_LNMTvsPT)) +       
  xlab("log2fc_LNMTvsPTlnm")+ylab("log2fc_LNMTvsPT")+      
  geom_point()+ geom_smooth(method="lm",colour='#764C29',fill='#E7E1D7') + 
  theme_bw()+      
  stat_cor(method = 'spearman', aes(x =log2fc_LNMTvsPTlnm, y =log2fc_LNMTvsPT))
ggMarginal(p2, type = "density", 
           xparams = list(fill = "#D3D3D3"),
           yparams = list(fill = "#D3D3D3"))


# Figure 3G
library(ggplot2)
library(ggstatsplot)
library(ggpubr)
kegg_c1 <- kegg_c1[1:5,] # top up gene result
ggplot(kegg_c1, aes(x=Description, y=-log10(p.adjust),
                    fill = "#FF6347")) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("-log10(FDR)") +
  theme(axis.text.x=element_text(hjust = 1,colour="black",
                                 size=10,angle = 90))
kegg_c2 <- kegg_c2[1:5,] # top down gene result
ggplot(kegg_c2, aes(x=Description, y=-log10(p.adjust),
                    fill = "#FF6347")) + 
  geom_bar(stat="identity",  position=position_dodge()) + 
  theme_bw() + ylab("-log10(FDR)") +
  theme(axis.text.x=element_text(hjust = 1,colour="black",
                                 size=10,angle = 90))



















