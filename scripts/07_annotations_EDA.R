##############################
## Data Annotations and EDA ##
##############################

##### import necessary libraries #####

library(Signac)
library(Seurat)
library(ggplot2)
library(GenomicRanges)
library(future)
library(reshape)
library(dplyr)
library(cowplot)
library(stringr)
library(hdf5r)
library(RColorBrewer)
library(glmGamPoi)
library(patchwork)
library(scCustomize)
library(Matrix)
library(dendextend)
library(cluster)
library(factoextra)
library(reshape2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)

#### Load Data ####

# Read RDS file storing pre-processed 2nd kidney PBMC scRNA-seq
sample <- readRDS("second_PBMC_batchcorr.rds")

#### Annotations ####

## Cell Types ##
Idents(sample) <- sample@meta.data$seurat_clusters
sample <- RenameIdents(sample, 
                       `0` = "Monocytes_CD14",
                       `1` = "Macrophage_1", 
                       `2` = "Macrophage_2",
                       `3` = "NK", 
                       `4` = "T_CD8", 
                       `5` = "T_CD4", 
                       `6` = "Naive_B", 
                       `7` = "Memory_B",
                       `8` = "T_Reg", 
                       `9` = "Megakaryocyte",
                       `10` = "NKT1",
                       `11` = "Megakaryocyte_Pro", 
                       `12` = "NKT2", 
                       `13` = "Erythrocytes_RBC",
                       `14` = "Monocytes_CD16",
                       `15` = "Plasma_B",
                       `16` = "Other")

Idents(sample) <- factor(Idents(sample),
                         levels = c(
                           "Monocytes_CD14",
                           "Monocytes_CD16",
                           "Macrophage_1",
                           "Macrophage_2",
                           "NK",
                           "NKT1",
                           "NKT2",
                           "T_CD8",
                           "T_CD4",
                           "T_Reg",
                           "Naive_B",
                           "Memory_B",
                           "Plasma_B",
                           "Megakaryocyte_Pro",
                           "Megakaryocyte",
                           "Erythrocytes_RBC",
                           "Other"))

cell_type <- sample@active.ident
names(cell_type) <- colnames(sample)
sample <- AddMetaData(object = sample, 
                      metadata = cell_type, 
                      col.name = "cell_type")

## Shortened Cell Type Names ##
Idents(sample) <- sample@meta.data$cell_type
sample <- RenameIdents(sample,
                       "Monocytes_CD14" = "mono-CD14",
                       "Monocytes_CD16" = "mono-CD16",
                       "Macrophage_1" = "M1",
                       "Macrophage_2" = "M2",
                       "Megakaryocyte" = "MGK",
                       "Megakaryocyte_Pro" = "MGK-P",
                       "Memory_B" = "M.B",
                       "Naive_B" = "N.B",
                       "Plasma_B" = "P.B",
                       "Erythrocytes_RBC" = "E-RBC",
                       "T_CD4" = "T-CD4",
                       "T_CD8" = "T-CD8",
                       "T_Reg" = "T-Reg")

Idents(sample) <- factor(Idents(sample),
                         levels = c("mono-CD14",
                                    "mono-CD16",
                                    "M1",
                                    "M2",
                                    "NK", "NKT1", "NKT2",
                                    "T-CD8","T-CD4",
                                    "T-Reg",
                                    "N.B", "M.B", "P.B",
                                    "MGK-P","MGK",
                                    "E-RBC", "Other"))

short <- sample@active.ident
names(short) <- colnames(sample)
sample <- AddMetaData(object = sample, metadata = short, col.name = "short")

## Shortened Time Annotations ##
Idents(sample) <- sample@meta.data$Condition
sample <- RenameIdents(sample, "PBMC_2nd_0h" = "0",
                       "PBMC_2nd_6h" = "6",
                       "PBMC_2nd_12h" = "12",
                       "PBMC_2nd_24h" = "24",
                       "PBMC_2nd_48h" = "48",
                       "PBMC_2nd_53h" = "53")

Idents(sample) <- factor(Idents(sample),
                         levels = c("0","6","12","24","48", "53"))

time <- sample@active.ident
names(time) <- colnames(sample)
sample <- AddMetaData(object = sample, metadata = time, col.name = "Time")

## Split Seurat Object by Cell Type ##
sample.list <- SplitObject(sample, split.by = "cell_type")

#### UMAP ####

## Figure 5.a ##
Idents(sample) <- sample@meta.data$short

fig_A <- 
  UMAPPlot(sample, label = TRUE, repel = TRUE) +
  theme(axis.title = element_text(size = 8, vjust = 0.5), 
        axis.text = element_text(size = 8)) + 
  NoLegend() 

#### Cell Type Proportions by Time ####

# Find percentage of cells from each cell type out of total cells at time point
cell_type_freq_by_condition <-  sample@meta.data %>%
  group_by(Time, short) %>%
  summarize(Frequency = n()) %>%
  mutate(Percentage = Frequency / sum(Frequency) * 100)

# 
cell_type_prop <- ggplot(cell_type_freq_by_condition, 
                         aes(x= Time, 
                             y = Percentage, 
                             fill = short, 
                             color = short)) + 
  geom_bar(stat = "identity", width = 0.5, color = "black")

## Figure 5.b ##
fig_B <- 
  cell_type_prop + 
  theme_classic() + 
  theme(axis.title = element_text(size = 10),
        axis.title.y = element_text(vjust = -1.5),
        axis.text = element_text(size = 10), 
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.margin = margin(0,0,0,0), 
        legend.box.margin = margin(c(-0.5, -0.5,-0.5,-0.5)), 
        legend.key.size = unit("0.25", "cm" )) +
  guides(color = guide_legend(override.aes = list(size = 1))) + 
  NoLegend()

#### QC Visualization ####

Idents(sample) <- sample@meta.data$Time
sample.time <- SplitObject(sample, split.by = "Time")
titles <- c("0h", "6h", "12h", "24h", "48h", "53h")
id <- 1
qc.plots <- c()

# Extract nFeature and nCount
for(i in sample.time){
  Idents(i) <- i@meta.data$short
  qc.plt1 <- VlnPlot(i, features = c("nFeature_RNA"),
                     split.by = "short",
                     flip = T, pt.size = 0) + 
    theme_classic() +
    ggtitle("nFeature_RNA") +
    theme(plot.title = element_text(size = 9, vjust = -1),
          axis.text.y = element_text(size = 9),
          axis.title = element_text(size = 9),
          axis.text.y.right = element_blank(), 
          axis.title.y.right = element_blank(), 
          axis.title.y.left = element_blank(), 
          axis.title.x.bottom = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank()) + NoLegend() +
    theme(plot.margin = margin(0.1,0.1,0.1,0.1, "mm"))
  
  qc.plt2 <- VlnPlot(i, features = c("nCount_RNA"),
                     split.by = "short",
                     flip = T, pt.size = 0) + 
    theme_classic() +
    ggtitle("nCount_RNA") +
    theme(plot.title = element_text(size = 9, vjust = -1),
          axis.text.x = element_text(size = 9, angle = 90),
          axis.text.y = element_text(size = 9),
          axis.title = element_text(size = 9),
          axis.text.y.right = element_blank(), 
          axis.title.y.right = element_blank(), 
          axis.title.y.left = element_blank(), 
          axis.title.x.bottom = element_blank()) + NoLegend() +
    theme(plot.margin = margin(0.1,0.1,0.1,0.1, "mm"))
  
  qc.plt <- ggarrange(qc.plt1, qc.plt2, 
                      ncol = 1, align = "v", heights = c(1,1.5))
  
  qc.plt <- annotate_figure(qc.plt, top = text_grob(titles[[id]], size = 9))
  qc.plots[[titles[[id]]]] <- qc.plt
  
  id <- id +1
  
}

## Supplementary Figure 5 ##
fig_S5 <- 
  ggarrange(plotlist = qc.plots, ncol = 2, nrow = 3, align = "hv", 
            labels = c("A", "B", "C", "D", "E", "F"))


#### Markers ####

## Supplementary Figure 6 ##
Idents(sample) <- sample@meta.data$short
fig_S6 <- 
  VlnPlot(sample, 
          features = c("LYZ", "CD14", "FCGR3A", "CD163", 
                       "GZMM",  "NKG7", "IL7R", "CD3D", "TRAC" ,"CCR7", "IL2RA",
                       "CD79A", "MS4A1", "TCL1A", "IGHA1", 
                       "GNG11", "PPBP", "HBA2", "HBB", 
                       "PRSS57", "STMN1", "SOX4"), 
          split.by = "short", stack = TRUE, flip = TRUE) + 
  NoLegend()

#### IFNG ####

store <- data.frame()

for(i in levels(cell_type)){
  
  cur <- DotPlot(sample.list[[i]], 
                 features = "IFNG", 
                 group.by = "Time", 
                 scale = FALSE)
  
  cur$data$features.plot <- as.factor(i)
  store <- rbind(store, cur$data)
}

sub <- rownames(sample)
sub <- sub[1:17]
frame <- DotPlot(sample, features = sub, group.by = "Time", scale = FALSE)
frame$data <- store

## Supplementary Figure 8 ##
fig_S8 <- 
  frame+ coord_flip() + 
  xlab("Cell Types") + 
  ylab("Hrs After Xenotransplantation") + 
  ggtitle("IFNG Expression")