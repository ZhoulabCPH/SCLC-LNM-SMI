###  Figure 1A
library(Seurat)
library(tidyverse)
library(ggrepel)
# TMA samples
length(unique(substr(PT$TMA,1,7)))
length(unique(substr(PT_LNM$TMA,1,7)))
length(unique(substr(LNMT$TMA,1,7)))



## Figure 1B
library(Nebulosa)
library(ggrastr)
library(Seurat)
library(tidyverse)
library(ggplot2)
dim_plot <- DimPlot(all_harmony8, reduction = "umap", 
                    group.by = "cellcomartment",raster=T,
                    cols = c("#E6C378","#79442D","#981868",
                             "#5F9EA0","#708090")) 
LabelClusters(dim_plot, id = "cellcomartment", repel = TRUE)
FeaturePlot(all_harmony8, min.cutoff = 0, max.cutoff = 1,
            features = c("PTPRC","CD68","JCHAIN","CD2", # Immune
                         "COL1A2","ACTA2",           # Fibroblast
                         "VWF","PECAM1",             # Endothelial
                         "EPCAM","KRT19" # Epithelial
            ), 
            reduction = "umap",
            pt.size = 0.5, cols = c("#F5F5F5","#8B0000"),
            alpha = c(0.8),raster = T) +
  theme(plot.title = element_text(hjust = 0.5))


# Figure 1C left
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggrepel)
dim_plot <- DimPlot(all_harmony8, reduction = "umap", 
                    group.by = "cell_ann", raster = TRUE,
                    cols = c("#E6C378","#3F51B5","#7B68EE","#9ACDFA",  
                             "#3498DB","#981868","#E26CD9","#708090")) + 
  theme(legend.position = "right") +
  theme_minimal() + 
  theme(axis.line = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_blank()) + 
  theme(legend.key = element_blank()) +  
  guides(color = guide_legend(override.aes = list(size = 3)))  
LabelClusters(dim_plot, id = "cell_ann", repel = TRUE)
# Figure 1C right
dim_plotPT <- DimPlot(PT, reduction = "umap", 
                          group.by = "cell_ann", raster = TRUE,
                          cols = c("#E6C378","#3F51B5","#7B68EE","#9ACDFA",  
                                   "#3498DB","#981868","#E26CD9","#708090")) + 
  theme(legend.position = "right") +
  theme_minimal() +  
  theme(axis.line = element_blank(),  
        panel.grid = element_blank(),  
        panel.border = element_blank()) + 
  theme(legend.key = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 3)))  
LabelClusters(dim_plotPT, id = "cell_ann", repel = TRUE)
dim_plotPT_LNM <- DimPlot(PT_LNM, reduction = "umap", 
                          group.by = "cell_ann", raster = TRUE,
                          cols = c("#E6C378","#3F51B5","#7B68EE","#9ACDFA",  
                                   "#3498DB","#981868","#E26CD9","#708090")) + 
  theme(legend.position = "right") +
  theme_minimal() +  
  theme(axis.line = element_blank(),  
        panel.grid = element_blank(), 
        panel.border = element_blank()) +  
  theme(legend.key = element_blank()) +  
  guides(color = guide_legend(override.aes = list(size = 3)))  
LabelClusters(dim_plotPT_LNM, id = "cell_ann", repel = TRUE)
dim_plotLNMT <- DimPlot(LNMT, reduction = "umap", 
                          group.by = "cell_ann", raster = TRUE,
                          cols = c("#E6C378","#3F51B5","#7B68EE","#9ACDFA",  
                                   "#3498DB","#981868","#E26CD9","#708090")) + 
  theme(legend.position = "right") +
  theme_minimal() +  
  theme(axis.line = element_blank(), 
        panel.grid = element_blank(),  
        panel.border = element_blank()) +  
  theme(legend.key = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 3)))  
LabelClusters(dim_plotLNMT, id = "cell_ann", repel = TRUE)
# Figure 1D
all_cell_num <- data.frame(num = c(dim(PT)[2],dim(PT_LNM)[2],
                                   dim(LNMT)[2]),
                           condition = c("PT","PT_LNM",
                                         "LNMT"))
PT_cell <- table(PT$cell_ann)/dim(PT)[2]
PT_LNM_cell <- table(PT_LNM$cell_ann)/dim(PT_LNM)[2]
LNMT_cell <- table(LNMT$cell_ann)/dim(LNMT)[2]
immune <- data.frame(per = as.numeric(c(PT_cell[2:5],PT_LNM_cell[2:5],
                                        LNMT_cell[2:5])),
                     cell = rep(c(names(PT_cell)[2:5]),3),
                     condition = c(rep("PT",4),rep("PT_LNM",4),
                                   rep("LNMT",4)))
malignant <- data.frame(per = as.numeric(c(PT_cell[c(1)],PT_LNM_cell[c(1)],
                                           LNMT_cell[c(1)])),
                        cell = rep(c(names(PT_cell)[c(1)]),3),
                        condition = c(rep("PT",1),rep("PT_LNM",1),
                                      rep("LNMT",1)))
other <- data.frame(per = as.numeric(c(PT_cell[c(6:8)],PT_LNM_cell[c(6:8)],
                                       LNMT_cell[c(6:8)])),
                    cell = rep(c(names(PT_cell)[c(6:8)]),3),
                    condition = c(rep("PT",3),rep("PT_LNM",3),
                                  rep("LNMT",3)))
ggplot(immune, aes(x = condition, y = per, fill = cell)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#3498DB","#7B68EE",
                               "#9ACDFA","#3F51B5"))+
  theme_bw() + ylim(0,0.3)+
  theme(axis.text.x = element_text(hjust = 1, colour = "black", size = 10, angle = 90), 
        axis.text.y = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10, colour = "black"), 
        legend.text = element_text(colour = "black", size = 10), 
        legend.title = element_text(colour = "black", size = 10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() ) +
  ylab("Percentage (%)")
ggplot(malignant, aes(x = condition, y = per, fill = cell)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#E5C378"))+
  theme_bw() + ylim(0,0.7)+
  theme(axis.text.x = element_text(hjust = 1, colour = "black", size = 10, angle = 90), 
        axis.text.y = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10, colour = "black"), 
        legend.text = element_text(colour = "black", size = 10), 
        legend.title = element_text(colour = "black", size = 10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() ) +
  ylab("Percentage (%)")
ggplot(other, aes(x = condition, y = per, fill = cell)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = c("#E26CD9",
                               "#981868","#708090"))+
  theme_bw() + ylim(0,0.3)+
  theme(axis.text.x = element_text(hjust = 1, colour = "black", size = 10, angle = 90), 
        axis.text.y = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10, colour = "black"), 
        legend.text = element_text(colour = "black", size = 10), 
        legend.title = element_text(colour = "black", size = 10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() ) +
  ylab("Percentage (%)")


# Figure 1E
library(Seurat)
library(tidyverse)
library(ggplot2)
all_harmony8$condition <- paste(all_harmony8$Patient_metastasis,
                                all_harmony8$Tissue_metastasis)
PT <- subset(all_harmony8, 
             subset = condition == "Primary Primary")
PT_LNM <- subset(all_harmony8, 
                 subset = condition == "Metastasis Primary")
LNMT <- subset(all_harmony8, 
               subset = condition == "Metastasis LNM")
# PT SMI  line 2 col 8 is our case
ImageDimPlot(PT,border.size = NA,fov = "YKKY0304_20240424_1",
             size = 0.05,group.by = "cell_ann",
             cols = c("#E6C378","#3F51B5","#7B68EE","#9ACDFA",
                      "#3498DB","#981868","#E26CD9","#708090"))
# PT_LNM SMI  line 6 col 1 is our case
ImageDimPlot(PT_LNM,border.size = NA,fov = "YKKY0304_20240424_2",
             size = 0.05,group.by = "cell_ann",
             cols = c("#E6C378","#3F51B5","#7B68EE","#9ACDFA",
                      "#3498DB","#981868","#E26CD9","#708090"))
# LNMT SMI  line 4 col 4 is our case
ImageDimPlot(LNMT,border.size = NA,fov = "YKKY0304_20240424_1",
             size = 0.05,group.by = "cell_ann",
             cols = c("#E6C378","#3F51B5","#7B68EE","#9ACDFA",
                      "#3498DB","#981868","#E26CD9","#708090"))

























