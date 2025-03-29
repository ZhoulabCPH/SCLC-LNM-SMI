
# Figure 5A and 5B
library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggvoronoi)
library(deldir)
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
###########
### PT type
###########
PT <- readRDS("PT_harmony8_afterQC_Doubcell_result.RDS")
PT_silhouette_plot_data <- readRDS("PT_silhouette.RDS")
PT_cell_metadata <- readRDS("PT_cell_metadata.RDS")
PT_cell_metadata$TMA <- PT$TMA[match(rownames(PT_cell_metadata),names(PT$cell_id))]
PT_cell_metadata$fov <- PT$fov[match(rownames(PT_cell_metadata),names(PT$cell_id))]
PT_OR_matrix <- readRDS("PT_OR_matrix.RDS")
PT_AdjPValue_matrix <- readRDS("PT_AdjPValue_matrix.RDS")
# Avg_Silhouette for dot plot
ggplot(PT_silhouette_plot_data, aes(x = Cluster, 
                                    y = Avg_Silhouette)) +
  geom_point() +
  geom_line() +
  xlab("Cluster") +
  ylab("Average Silhouette Width") +
  ggtitle("") +
  scale_x_continuous(breaks = seq(2, 15, 1)) +  # 
  theme_bw()
# Voronoi plot for case 1
voronoi_data <- PT_cell_metadata[which(PT_cell_metadata$TMA == "924838-7"),]
voronoi_data <- voronoi_data[which(voronoi_data$fov == "Slide1_47"),]
x_min <- min(voronoi_data$x_coord)
x_max <- max(voronoi_data$x_coord)
y_min <- min(voronoi_data$y_coord)
y_max <- max(voronoi_data$y_coord)
boundary <- data.frame(
  x = c(x_min, x_max, x_max, x_min),
  y = c(y_min, y_min, y_max, y_max))
ggplot(voronoi_data, aes(x = x_coord, y = y_coord, 
                         fill = as.factor(cellular_neighborhood))) +
  geom_voronoi(aes(fill = as.factor(cellular_neighborhood)), 
               outline = boundary, color = NA) +
  scale_fill_manual(values = c("1" = "#9E9E9A","5" = "#E9C46A",
                               "7" = "#FF9AA2","3" = "#A0E7E5",
                               "9" = "#7B68EE","11" = "#0D98BA")) +
  theme_void() + coord_fixed()
# Heatmap for 5B
PT_log_OR_matrix <- log2(PT_OR_matrix+0.5) # normalized
PT_AdjPValue_matrix <- PT_AdjPValue_matrix[,colnames(PT_log_OR_matrix)]
PT_log_OR_matrix[PT_AdjPValue_matrix > 0.05] <- 0  # normalized
PT_log_OR_matrix <- PT_log_OR_matrix[,name]
pheatmap(PT_log_OR_matrix,fontsize=6,
         cellwidth = 13,cellheight = 13,
         color  = colorRampPalette(c("#7F3C8D","#F2F2F2","#FFCB05"))(100),
         border_color = "grey60",
         breaks = seq(-2, 2, length.out = 101),
         scale = "row",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)
###############
### PT-LNM type
###############
PT_LNM <- readRDS("PT_LNM_harmony8_afterQC_Doubcell_result.RDS")
PT_LNM_silhouette_plot_data <- readRDS("PT_LNM_silhouette.RDS")
PT_LNM_cell_metadata <- readRDS("PT_LNM_cell_metadata.RDS")
PT_LNM_cell_metadata$TMA <- PT_LNM$TMA[match(rownames(PT_LNM_cell_metadata),names(PT_LNM$cell_id))]
PT_LNM_cell_metadata$fov <- PT_LNM$fov[match(rownames(PT_LNM_cell_metadata),names(PT_LNM$cell_id))]
PT_LNM_OR_matrix <- readRDS("PT_LNM_OR_matrix.RDS")
PT_LNM_AdjPValue_matrix <- readRDS("PT_LNM_AdjPValue_matrix.RDS")
# Avg_Silhouette for dot plot
ggplot(PT_LNM_silhouette_plot_data, aes(x = Cluster, 
                                        y = Avg_Silhouette)) +
  geom_point() +
  geom_line() +
  xlab("Clusters") +
  ylab("Average Silhouette Width") +
  ggtitle("") +
  scale_x_continuous(breaks = seq(2, 15, 1)) + 
  theme_bw()
# Voronoi plot for case 2
voronoi_data <- PT_LNM_cell_metadata[which(PT_LNM_cell_metadata$TMA == "682067-6"),]
voronoi_data <- voronoi_data[which(voronoi_data$fov == "Slide2_173"),]
x_min <- min(voronoi_data$x_coord)
x_max <- max(voronoi_data$x_coord)
y_min <- min(voronoi_data$y_coord)
y_max <- max(voronoi_data$y_coord)
boundary <- data.frame(
  x = c(x_min, x_max, x_max, x_min),
  y = c(y_min, y_min, y_max, y_max))
ggplot(voronoi_data, aes(x = x_coord, y = y_coord, 
                         fill = as.factor(cellular_neighborhood))) +
  geom_voronoi(aes(fill = as.factor(cellular_neighborhood)), 
               outline = boundary, color = NA) +
  scale_fill_manual(values = c("1" = "#FF9AA2","3" = "#9E9E9A",
                               "5" = "#E9C46A","10" = "#9E9E9A",
                               "11" = "#FF9AA2","6" = "#7B68EE",
                               "9" = "#0D98BA")) +
  theme_void() + coord_fixed()
# Heatmap for 5B
PT_LNM_log_OR_matrix <- log2(PT_LNM_OR_matrix+0.5)
PT_LNM_AdjPValue_matrix <- PT_LNM_AdjPValue_matrix[,colnames(PT_LNM_log_OR_matrix)]
PT_LNM_log_OR_matrix[PT_LNM_AdjPValue_matrix > 0.05] <- 0
PT_LNM_log_OR_matrix <- PT_LNM_log_OR_matrix[,name]
pheatmap(PT_LNM_log_OR_matrix,fontsize=6,
         cellwidth = 13,cellheight = 13,
         color  = colorRampPalette(c("#7F3C8D","#F2F2F2","#FFCB05"))(100),
         border_color = "grey60",
         breaks = seq(-2, 2, length.out = 101),
         scale = "row",
         cluster_cols = F, cluster_rows = F,
         show_rownames = T, show_colnames = T)


# Figure 5C
library(circlize)
library(ggplot2)
library(ggforce)
library(ggpubr)
PT_rate <- readRDS("PT_CNrate.RDS")
PT_LNM_rate <- readRDS("PT_LNM_CNrate.RDS")
# PT
PT_percentage <- c()
for (i in 1:11){
  weizhi <- which(PT_rate$cellular_neighborhood == i)
  PT_percentage <- c(PT_percentage,
                     sum(PT_rate$percentage[weizhi])/sum(PT_rate$percentage))
}
PT_mat <- data.frame(percentage = PT_percentage,
                     cn = as.character(c(1:11)),
                     type = c("Undefined","Tumour_compartment",
                              "Immune hotspot-2","Tumour_compartment",
                              "Tumour_compartment","Exh_tumour_compartment",
                              "Exh_tumour_compartment","Exh_tumour_compartment",
                              "Immune hotspot-1","Tumour_compartment","Vascular"))
PT_mat <- PT_mat[order(PT_mat$type),]
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+#添加颜色
  scale_fill_manual(values = c('#FF9AA2','#7B68EE','#A0E7E5',
                               '#0D98BA','#E9C46A','#9E9E9A'))+
  geom_arc_bar(data=PT_mat,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=percentage,fill=type)
  )
# PT-LNM
PT_LNM_percentage <- c()
for (i in 1:13){
  weizhi <- which(PT_LNM_rate$cellular_neighborhood == i)
  PT_LNM_percentage <- c(PT_LNM_percentage,
                         sum(PT_LNM_rate$percentage[weizhi])/sum(PT_LNM_rate$percentage))
}
PT_LNM_mat <- data.frame(percentage = PT_LNM_percentage,
                         cn = as.character(c(1:13)),
                         type = c("Exh_tumour_compartment","Tumour_compartment",
                                  "Undefined","Tumour_compartment",
                                  "Tumour_compartment","Immune hotspot-1",
                                  "Tumour_compartment","Tumour_compartment",
                                  "Vascular","Undefined",
                                  "Exh_tumour_compartment","Exh_tumour_compartment",
                                  "Tumour_compartment"))
PT_LNM_mat <- PT_LNM_mat[order(PT_LNM_mat$type),]
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+#添加颜色
  scale_fill_manual(values = c('#FF9AA2','#7B68EE','#0D98BA',
                               '#E9C46A','#9E9E9A'))+
  geom_arc_bar(data=PT_LNM_mat,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=percentage,fill=type)
  )
options(warn = -1)
clinical <- read.csv("CNs_percentage_Patient_clinical_information.csv",
                     stringsAsFactors = F,
                     row.names = 1,check.names = F)
clinical_box <- data.frame(rate = c(clinical$tumor,clinical$Exh_tumo,
                                    clinical$immhos1,clinical$immhos2,
                                    clinical$bcell),
                           CN = c(rep("Tumour compartment",nrow(clinical)),
                                  rep("Exh. tumour compartment",nrow(clinical)),
                                  rep("Immune hotspot-1",nrow(clinical)),
                                  rep("Immune hotspot-2",nrow(clinical)),
                                  rep("B ell enriched",nrow(clinical))),
                           meta = rep(clinical$Patient_metastasis,5))
ggplot() +
  geom_half_boxplot(data = clinical_box, 
                    aes(x=CN, y=rate, fill = meta),
                    width = 0.3, position = position_dodge(0.7), 
                    alpha=1, outlier.shape = NA)+
  geom_half_point(data = clinical_box, 
                  aes(x=CN, y=rate, color=meta), 
                  shape = 16, alpha=1, size = 1,
                  transformation = position_jitter(width = 0.1,height = 0))+
  scale_fill_manual(values = c('#d27f69','#97b6a0'))+
  scale_color_manual(values = c('#d27f69','#97b6a0'))+
  stat_compare_means(data = clinical_box,
                     aes(x=CN, y=rate, fill=meta),
                     method = "wilcox.test")+   theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.position="right") + labs(x='', y="rate")


# Figure 5D
library(survival)
library(survminer)
library(ggplot2)
clinical <- read.csv("CNs_percentage_Patient_clinical_information.csv",
                     stringsAsFactors = F,
                     row.names = 1,check.names = F)
clinical$immhos1_cluster_new <- clinical$immhos1
clinical$immhos1_cluster_new[clinical$immhos1 > median(clinical$immhos1)] <- "high"
clinical$immhos1_cluster_new[clinical$immhos1 <= median(clinical$immhos1)] <- "low"
clinical$immhos1_cluster_new[clinical$immhos1 == 0] <- "None"
# immhos1_cluster_new
# OS
surv <- survfit(Surv(OS_time,OS_status)~immhos1_cluster_new,data = clinical)
survdiff(Surv(OS_time,OS_status)~immhos1_cluster_new,data = clinical)
pairwise_survdiff(Surv(OS_time, OS_status) ~ immhos1_cluster_new, 
                  p.adjust.method = "none",data = clinical)
summary(surv)
summary(coxph(Surv(OS_time,OS_status)~immhos1_cluster_new,data = clinical))
plot_os <- ggsurvplot(surv,pval = TRUE,
                      risk.table = TRUE, # Add risk table
                      risk.table.col = "strata", # Change risk table color by groups
                      palette = c("#FFB677","#98C6E1","#0A1931"),
                      legend = c(2,0.5), 
                      legend.title = "", 
                      conf.int = F,
                      legend.labs = c("High","low","None"), # 指定图例分组标签
                      xlab = "Time (months)",
                      ylab = "Overall Survival")
plot_os$plot <- plot_os$plot +
  geom_vline(xintercept = 36, linetype = "dashed", color = "grey", size = 0.7) +  
  geom_vline(xintercept = 60, linetype = "dashed", color = "grey", size = 0.7)  
plot_os
# DFS
surv <- survfit(Surv(DFS_time,DFS_status)~immhos1_cluster_new,data = clinical)
survdiff(Surv(DFS_time,DFS_status)~immhos1_cluster_new,data = clinical)
pairwise_survdiff(Surv(DFS_time, DFS_status) ~ immhos1_cluster_new, 
                  p.adjust.method = "none",data = clinical)
summary(surv)
summary(coxph(Surv(DFS_time,DFS_status)~immhos1_cluster_new,data = clinical))
plot_dfs <- ggsurvplot(surv,pval = TRUE,
                       risk.table = TRUE, # Add risk table
                       risk.table.col = "strata", # Change risk table color by groups
                       palette = c("#FFB677","#98C6E1","#0A1931"),
                       legend = c(2,0.5), 
                       legend.title = "", 
                       conf.int = F,
                       legend.labs = c("High","low","None"), 
                       xlab = "Time (months)",
                       ylab = "DFS")
plot_dfs$plot <- plot_dfs$plot +
  geom_vline(xintercept = 36, linetype = "dashed", color = "grey", size = 0.7) +  
  geom_vline(xintercept = 60, linetype = "dashed", color = "grey", size = 0.7)  
plot_dfs


# Figure 5E
library(forestplot)
data2 <- read.csv("CNs_os_多cox.csv",stringsAsFactors = F)
forestplot(as.matrix(data2[,c(1,4,9)]),data2$exp.coef.,data2$lower..95,
           data2$upper..95,zero = 1,xlog = F,
           clip = c(0,5),
           colgap = unit(6,"mm"),graphwidth=unit(70,"mm"),
           lineheight = unit(0.3,"cm"),graph.pos = 3,
           col = fpColors(box="black", lines="black", zero = "black"),boxsize = 0.2,
           ci.vertices = T,ci.vertices.height = 0.2,
           lty.ci = 7,lwd.zero=0.4,lwd.ci = 3,
           txt_gp=fpTxtGp(label = gpar(cex=0.8),
                          ticks = gpar(cex=0.8)) ,
           is.summary=c(TRUE,rep(FALSE,100)))
data2 <- read.csv("CNs_dfs_多cox.csv",stringsAsFactors = F)
forestplot(as.matrix(data2[,c(1,4,9)]),data2$exp.coef.,data2$lower..95,
           data2$upper..95,zero = 1,xlog = F,
           clip = c(0,4),
           colgap = unit(6,"mm"),graphwidth=unit(70,"mm"),
           lineheight = unit(0.3,"cm"),graph.pos = 3,
           col = fpColors(box="black", lines="black", zero = "black"),boxsize = 0.2,
           ci.vertices = T,ci.vertices.height = 0.2,
           lty.ci = 7,lwd.zero=0.4,lwd.ci = 3,
           txt_gp=fpTxtGp(label = gpar(cex=0.8),
                          ticks = gpar(cex=0.8)) ,
           is.summary=c(TRUE,rep(FALSE,100)))











