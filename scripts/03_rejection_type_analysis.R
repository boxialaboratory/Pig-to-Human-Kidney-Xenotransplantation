#################################
##### Types of Rejection ########
#################################
library(tidyverse)
library(reshape2)
library(RColorBrewer)

##### Figure 3 #####
##### Bulk Data for Expression over Time #####
bulk <- read.csv("data/2nd_kidney_bulk.csv")
colnames(bulk) <- gsub(".*T|\\.counts\\.rev\\.txt", "", colnames(bulk))
rownames(bulk) <- bulk$gene
bulk <- bulk[,-c(1)]

# Calculate average of T0.1 and T0.2
bulk$T0 <- rowMeans(bulk[, c("0.1", "0.2")])
bulk <- bulk[, setdiff(names(bulk), c("0.1", "0.2"))]

# Calculate log10 fold change across different time points for each gene.
# Here, we're calculating the fold change relative to the average T0 value we just calculated.
bulk <- log2(bulk / bulk[,"T0"])
bulk$gene <- rownames(bulk)
bulk <- melt(bulk, id.vars="gene")
bulk$variable <- factor(bulk$variable, levels = c("T0", "12", "24", "48"))
# Antibody-mediated rejection genes
genes_abmr <- c(
	"CXCL10", "CXCL11", "PLA1A", "CCL4d", "CD160", "YME1L1", "FGFBP2", "GNLY", "CX3CR1", 
	"KLRD1", "KLRF1", "SH2D1B", "TRDV3", "CDH5", "CDH13", "COL13A1", "DARC", "ECSCR", 
	"GNG11", "ICAM2", "MALL", "PECAM1", "PGM5", "RAMP3", "RAPGEF5", "ROBO4", "TM4SF18", 
	"VWF", "KLF4", "PPM1F"
)
# T-cell mediated rejection genes
genes_tcmr <- c(
	"ANKRD22", "APOL2", "GBP5", "RARRES3", "LAP3", "CTLA4", "ICOS", "CD96", "IFNG", "LAG3", 
	"SIRPG", "BTLA", "LCP2", "DUSP2", "MYB", "CD8A", "IL12RB1", "AOAH", "TNFSF8", "CD86", 
	"CD274", "ADAMDEC1", "CD84", "FGD2", "IL21R", "LAIR1", "PHEX", "SLA", "SLAMF8", "CD72"
)

# Define the order of importance for genes of each category (Replace with your own ordering)
gene_order_abmr <- c("CXCL10", "CXCL11", "PLA1A", "CCL4d", "CD160", "YME1L1", "FGFBP2", "GNLY", "CX3CR1", 
							"KLRD1", "KLRF1", "SH2D1B", "TRDV3", "CDH5", "CDH13", "COL13A1", "DARC", "ECSCR", 
							"GNG11", "ICAM2", "MALL", "PECAM1", "PGM5", "RAMP3", "RAPGEF5", "ROBO4", "TM4SF18", 
							"VWF", "KLF4", "PPM1F")
gene_order_tcmr <- c("ANKRD22", "APOL2", "GBP5", "RARRES3", "LAP3", "CTLA4", "ICOS", "CD96", "IFNG", "LAG3", 
							"SIRPG", "BTLA", "LCP2", "DUSP2", "MYB", "CD8A", "IL12RB1", "AOAH", "TNFSF8", "CD86", 
							"CD274", "ADAMDEC1", "CD84", "FGD2", "IL21R", "LAIR1", "PHEX", "SLA", "SLAMF8", "CD72")

check_genes_and_assign_colors <- function(genes, color, category, gene_order){
	present_genes <- genes[genes %in% bulk$gene] # Check which genes are present in the dataset
	present_genes <- present_genes[match(present_genes, gene_order)] # Sort genes according to importance
	color_palette <- rev(colorRampPalette(brewer.pal(9, color))(length(present_genes))) # Create and reverse color palette
	
	# Filter for selected genes and assign colors
	subset_data <- bulk[bulk$gene %in% present_genes,]
	subset_data$gene <- factor(subset_data$gene, levels = present_genes)
	subset_data$color <- color_palette[match(subset_data$gene, present_genes)]
	subset_data$category <- category
	
	# Return a list containing the subset data and the gene-color mapping
	return(list("data" = subset_data, "gene_color_mapping" = setNames(color_palette, present_genes)))
}

# Get the data and the gene-color mapping for each type of rejection
abmr <- check_genes_and_assign_colors(genes_abmr, "Reds", "Antibody-mediated Rejection", gene_order_abmr)
tcmr <- check_genes_and_assign_colors(genes_tcmr, "Blues", "T-Cell Mediated Rejection", gene_order_tcmr)

# define function for plotting all rejection-related markers
plotBulkCategories <- function(data, gene_color_mapping, title){
  ggplot(data, aes(x = as.numeric(variable), y = value, group = gene, color = gene)) +
    geom_line(size = 1.5) +
    theme_classic() +
    
    ylim(-2,6) +
    scale_color_manual(values = gene_color_mapping) +
    theme(legend.position = "bottom",
          legend.text = element_text(size=12),
          legend.title = element_blank()) +
    labs(x = "Hours after transplantation", y = "log2FC gene expression", color = "Gene") +
    ggtitle(title)
}

# Figure 3 A
plotBulkCategories(abmr$data, 
						 abmr$gene_color_mapping, 
						 "Antibody-mediated Rejection Over Time")
# Figure 3 B
plotBulkCategories(tcmr$data, 
						 tcmr$gene_color_mapping, 
						 "T-Cell Mediated Rejection Over Time")
##### Endothelial Activation ####
# Figure 3 C
VlnPlot(subset(kidneys,
subset = cell_types == "Endo"),
features = c(
	"CDH5",
	"CDH13", 
	"PECAM1",
	"RAMP3", 
	"TM4SF18"),
group.by = "Condition",
split.by = "cell_types",
pt.size = 0.05,
combine = TRUE,
flip = TRUE,
same.y.lims = TRUE,
stack = TRUE) + 
	scale_fill_manual(values = "pink") +
	theme_classic()+
	theme(axis.title = element_blank())+
	stat_compare_means(aes(group = Condition), 
label = "p.signif",
method = "wilcox.test", 
comparisons = list(c("kidney_1st_untreated_rt","kidney_1st_treated"), c("kidney_2nd_untreated_rt","kidney_2nd_treated"))) +
	ylim(0,6)

# Figure 3 D
VlnPlot(Immune1,
		  features = c("ANKRD22",
		  		"SLAMF8"),
		  group.by = "Immune_Idents",
		  split.by = "Condition",
		  pt.size = 0.01,
		  stack = TRUE,
		  combine = TRUE,
		  flip = TRUE,
		  same.y.lims = TRUE,
		  cols = c("#FDD6D5","#FF7976"))+
	ylim(NA, 3)+
	theme_classic2() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1),
			axis.title = element_blank())+
	stat_compare_means(aes(group = Immune_Idents), 
							 label = "p.signif", 
							 method = "wilcox.test", 
							 comparisons = list(c("Porcine_NKT","Human_NK"),
							 						 c("Porcine_Macrophage","Human_Macrophage")))

VlnPlot(Immune2,
		  features = c("ANKRD22",
		  				 "SLAMF8"),
		  group.by = "Immune_Idents",
		  split.by = "Condition",
		  pt.size = 0.01,
		  stack = TRUE,
		  combine = TRUE,
		  flip = TRUE,
		  same.y.lims = TRUE,
		  cols = c("#8FD0FF","#5081ee"))+
	ylim(NA, 3)+
	theme_classic2() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1),
			axis.title = element_blank())+
	stat_compare_means(aes(group = Immune_Idents), label = "p.signif", method = "wilcox.test", 
							 comparisons = list(c("Porcine_NKT","Human_NK"),
							 						 c("Porcine_Macrophage","Human_Macrophage")))


