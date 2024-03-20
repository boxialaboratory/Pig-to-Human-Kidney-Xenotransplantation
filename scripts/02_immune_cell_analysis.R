###########################
## Immune Cell Analysis ##
###########################

##### import necessary libraries #####
library(ggplot2)
library(reshape2)
library(dplyr)
library(readr)
library(ggpubr)
library(xlsx)
library(Seurat)
library(patchwork)
library(RColorBrewer)
library(gridExtra)
##### cell type identification #####
# UMAP Generation
DimPlot(kidneys,
group.by = "origin",
cols = c("#E87D72","#999999"))+ 
theme(plot.title = element_blank()) 

# since most human cells are immune cells, zoom out to all immune cells
Immune <- subset(kidneys, subset = cell_types == "Immune")
# Subset into clusters
DefaultAssay(Immune) <- "RNA"
Immune <- NormalizeData(Immune, normalization.method = "LogNormalize", scale.factor = 10000)
Immune <- FindVariableFeatures(Immune, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(Immune)
Immune <- ScaleData(Immune, features = all.genes)
Immune <- RunPCA(Immune, features = VariableFeatures(object = Immune))
Immune <- FindNeighbors(Immune, dims = 1:10)
Immune <- FindClusters(Immune, resolution = 0.1)
Immune <- RunUMAP(Immune, dims = 1:10)

DefaultAssay(Immune) <- "SCT"
DimPlot(Immune,
reduction = "umap",
label = T,
split.by = "origin",
ncol = 2)

# Figure 2 A
DimPlot(Immune,
split.by = "Condition",
group.by = "origin",
ncol = 2) + theme(plot.title = element_blank(),
strip.text = element_blank(),
axis.title = element_blank()) + NoLegend()
DimPlot(Immune,
group.by = "origin") + theme(plot.title = element_blank())

# Supplementary Figure 4 A
FeaturePlot(Immune,
feature = c("CD163", "TYROBP"),
label = T,
split.by = "origin",
ncol = 3)
# Supplementary Figure 4 B
FeaturePlot(Immune,
feature = c("NKG7", "GNLY"),
label = T,
split.by = "origin")
# Supplementary Figure 4 C
FeaturePlot(Immune,
feature = c("CD3E", "CD3G"),
label = T,
split.by = "origin")

# Based on expression of marker genes, rename identities of the clusters
# Subset into Human Macrophages, Human NK Cells, Porcine T cells
Immune_Idents <- character(ncol(Immune))
# Iterate over each cell and assign the appropriate value to new_idents
for (i in 1:ncol(Immune)) {
	if (Immune[["seurat_clusters"]][i,] %in% c(0,3) && Immune[["origin"]][i,] == "human") {
		Immune_Idents[i] <- "Human_NK"
	} else if (Immune[["seurat_clusters"]][i,] %in% c(2) && Immune[["origin"]][i,] == "human") {
		Immune_Idents[i] <- "Human_Macrophage"
	} else if (Immune[["seurat_clusters"]][i,] %in% c(0) && Immune[["origin"]][i,] == "porcine") {
		Immune_Idents[i] <- "Porcine_NKT"
	} else if (Immune[["seurat_clusters"]][i,] %in% c(1) && Immune[["origin"]][i,] == "porcine") {
		Immune_Idents[i] <- "Porcine_T"
	} else if (Immune[["seurat_clusters"]][i,] %in% c(4) && Immune[["origin"]][i,] == "porcine") {
		Immune_Idents[i] <- "Porcine_Macrophage"
	}
}

Immune$Immune_Idents <- Immune_Idents
Immune <- SetIdent(Immune, value = "Immune_Idents")

# Figure 2 B
DimPlot(Immune) + theme(axis.title = element_blank())

##### Human Cell Infiltration %%%%%
# Read the csv data
bulk_ratio <- read.xlsx("data/human_to_pig_readsratio_wo_duplicate.xlsx", sheetIndex = 1)
markers <- c(# macrophage markers
"LYZ","CD163", 
"CCR5", "TYROBP",
# NK cell markers
"NKG7","GNLY","KLRD1",
"KLRE1","KLRF1","IFNG")

ifng_inducible <- c("IFNG", #interferon-gamma
 "CXCL9","CXCL10","CXCL11")#IFNG-induced cytokines)

bulk_markers <- bulk_ratio %>%
filter(hp_ortho_genes %in% markers) %>%
mutate(time_0h = (time_0h1 + time_0h2) / 2)

bulk_ifng_inducible <- bulk_ratio %>%
filter(hp_ortho_genes %in% ifng_inducible) %>%
mutate(time_0h = (time_0h1 + time_0h2) / 2)

#  Melt the dataframe to make it suitable for ggplot2, keep the genes and the times
bulk_markers_melt <- melt(bulk_markers, id.vars = "hp_ortho_genes", 
measure.vars = c("time_0h", "time_12h", "time_24h", "time_48h"))

bulk_ifng_inducible_melt <- melt(bulk_ifng_inducible, id.vars = "hp_ortho_genes", 
measure.vars = c("time_0h", "time_12h", "time_24h", "time_48h"))

# Figure 2 C
ggplot(bulk_markers_melt,  aes(x = variable, y = factor(hp_ortho_genes, 
levels = rev(markers)), 
fill = value)) +
geom_tile() +
geom_text(aes(label = round(value, 3)), color = "black", size = 5) +
scale_fill_gradientn(colors = c("lightblue","#FAE06E","#FA706A"),
limits=c(0, 2),
oob = scales::squish) +
labs(x = "Time", y = "Gene", fill = "Reads Ratio") +
theme_classic()+
theme(legend.position = "bottom")+ 
theme(legend.title = element_blank(),
axis.title = element_blank())
# Figure 2 F
ggplot(bulk_ifng_inducible_melt, aes(x = variable, y = factor(hp_ortho_genes, 
levels = rev(ifng_inducible)), 
fill = value)) +
geom_tile() +
geom_text(aes(label = round(value, 3)), color = "black", size = 5) +
scale_fill_gradientn(colors = c("lightblue","#FAE06E","#FA706A"),
limits=c(0, 2),
oob = scales::squish) +
labs(x = "Time", y = "Gene", fill = "Reads Ratio") +
theme_classic()+
theme(legend.position = "bottom") + 
theme(legend.title = element_blank(),
axis.title = element_blank())

Immune1 <- subset(Immune, subset = Condition %in% c("kidney_1st_untreated_rt","kidney_1st_treated"))
Immune2 <- subset(Immune, subset = Condition %in% c("kidney_2nd_untreated_rt","kidney_2nd_treated"))
# Specify the order of cell types
cell_type_order <- c("Porcine_T", "Porcine_NKT", "Human_NK", "Porcine_Macrophage", "Human_Macrophage")
# Reorder the levels of the cell type factor variable in your Seurat object
Immune1$Immune_Idents <- factor(Immune1$Immune_Idents, levels = cell_type_order)
Immune2$Immune_Idents <- factor(Immune2$Immune_Idents, levels = cell_type_order)

# Figure 2 D
VlnPlot(Immune1,
 features = ifng_inducible,
 group.by = "Immune_Idents",
 split.by = "Condition",
 pt.size = 0.01,
 stack = TRUE,
 combine = TRUE,
 flip = TRUE,
 same.y.lims = TRUE,
 cols = c("#FDD6D5","#FF7976","#8FD0FF","#5081ee"))+
 ylim(NA, 8)+
 theme_classic2() +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1),
 axis.title = element_blank())+
stat_compare_means(aes(group = Immune_Idents), label = "p.signif", method = "wilcox.test", 
comparisons = list(c("Porcine_NKT","Human_NK"), c("Porcine_Macrophage","Human_Macrophage")))

VlnPlot(Immune2,
 features = ifng_inducible,
 group.by = "Immune_Idents",
 split.by = "Condition",
 pt.size = 0.01,
 stack = TRUE,
 combine = TRUE,
 flip = TRUE,
 same.y.lims = TRUE,
 cols = c("#FDD6D5","#FF7976","#8FD0FF","#5081ee"))+
 ylim(NA, 7)+
 theme_classic2() +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1), axis.title = element_blank())+
 stat_compare_means(aes(group = Immune_Idents), label = "p.signif", method = "wilcox.test", 
comparisons = list(c("Porcine_NKT","Human_NK"),
c("Porcine_Macrophage","Human_Macrophage")))

##### Longitudinal bulk-RNA-seq #####
library(tidyverse)
library(reshape2)
bulk <- read.csv("data/2nd_kidney_bulk.csv")
# Preprocess column names to extract time points
# Now we're distinguishing between the T0.1 and T0.2 replicates
colnames(bulk) <- gsub(".*T|\\.counts\\.rev\\.txt", "", colnames(bulk))
rownames(bulk) <- bulk$gene
bulk <- bulk[,-c(1)]

# Calculate average of T0.1 and T0.2
bulk$T0 <- rowMeans(bulk[, c("0.1", "0.2")])
bulk <- bulk[, setdiff(names(bulk), c("0.1", "0.2"))]

# Here, we're calculating the fold change relative to the average T0 value we just calculated.
bulk <- log2(bulk / bulk[,"T0"])
bulk$gene <- rownames(bulk)
bulk <- melt(bulk, id.vars="gene")
bulk$variable <- factor(bulk$variable, levels = c("T0", "12", "24", "48"))
plotBulk <- function(genes){
# Filter for selected genes
subset_data <- bulk[bulk$gene %in% genes,]
subset_data$gene <- factor(subset_data$gene, levels = genes)
# Plot log2FC over time
ggplot(subset_data, aes(x = as.numeric(variable), y = value, group = gene, color = gene)) +
geom_line(size = 1.5) +
labs(x = "Hours after xenotransplantation", y = "log2FC Gene Expression") +
theme_classic() +
theme(legend.position="bottom",
legend.text = element_text(size=12),
legend.title = element_blank())
}

# Figure 2 E
plotBulk(ifng_inducible)+
scale_color_manual(values = c("#FDAE6B","#238B45","#41AB5D", "#A1D99B"))+
theme(plot.title = element_blank())
