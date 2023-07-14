###########################################
## Co-Expressed Gene Clustering Analysis ##
###########################################

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


####Percentage Expression Analysis####

#Run across all cell types
#Extract percentage expression for all genes in each cell types
#Subset genes with >20% expression
#Remove MT and RP genes
gene_list.list <- c()
pct.expression.list <- c()
genes.20.list <- c()
genes.list <- c()
cell_types <- c(names(sample.list))

id <- 1
for (i in sample.list){
  
  pct.expression <- Percent_Expressing(i, 
                                       features = rownames(i), 
                                       entire_object = TRUE)
  
  pct.expression.list[[id]] <- pct.expression
  genes_20 <- pct.expression %>% filter_at(1, all_vars(.>20))
  genes.20.list[[id]] <- genes_20
  genes <- tibble::rownames_to_column(genes_20, "gene")
  MT_RB.list <- genes[grepl("^MT|^RP", genes$gene),]
  genes <- genes %>% subset(!(gene %in% MT_RB.list$gene))
  genes.list[[id]]<- genes
  
  id <- id +1
  
}

# Name gene list extracted for each cell types by cell type name 
names(genes.list) <- cell_types

#### Gene Clustering ####
cluster.info.list <- c()
cluster.means.list <- c()

for(i in names(sample.list)){
  
  cur_cell_type <- sample.list[[i]]
  cur_genes_20 <- genes.list[[i]]
  
  # Extract normalized count matrix from SCT assay
  mtx <- GetAssayData(cur_cell_type, assay = "SCT", slot = "count")
  
  # Store normalized count matrix in new Seurat object
  data <- CreateSeuratObject(counts = mtx, assay = "SCT")
  
  # Add meta-data for time point annotations
  data <- AddMetaData(data, cur_cell_type@meta.data$Time, col.name = "Time")
  data@meta.data$Time <- factor( data@meta.data$Time,
                                 levels = c("0", "6", "12", "24", "48", "53"))
  
  DefaultAssay(data) <- "SCT"
  
  # Normalize by gene (feature) using CLR, make genes comparable for clustering
  data <- NormalizeData(data, normalization.method = "CLR", margin = 1)
  
  # Scale data
  data <- ScaleData(data, assay = "SCT", features = rownames(data))
  
  # Find variable features
  Idents(data) <- data@meta.data$Time
  data <- FindVariableFeatures(data, nfeatures = 2000)
  
  # Extract variable features that are expressed by 20% of cells
  genes <- VariableFeatures(data)
  genes <- genes %>% subset(genes %in% cur_genes_20$gene)
  
  # Find average expression of gene per time point
  xpr <- AverageExpression(data, assay = "SCT", slot = "scale.data", features = genes, group.by = "Time")
  xpr <- as.data.frame(xpr$SCT)
  
  # Generate correlation and distance matrix for genes
  cur.df <- xpr
  cor.df <- as.data.frame(t(cur.df))
  cor.df <- tibble::rownames_to_column(cor.df, "time")
  cor_mat <- dplyr::select(cor.df, -time) %>% cor(use="pairwise.complete.obs")
  dist_mat <- as.dist(1-cor_mat)
  
  # Create dendrogram by hierarchical clustering of genes
  cur.df<- xpr
  tree <- hclust(dist_mat, method = "complete")
  dend_plot <- plot(tree, cex = 0.1)
  dend <- as.dendrogram(tree)
  
  # Generate gene clusters by cutting dendrogram at height 1.5
  hclusters <- cutree(dend, h = 1.5)
  
  # Render cluster information
  hclusters.df <- data.frame(gene = names(hclusters), cluster = as.factor(hclusters))
  
  # Annotate cluster information onto scaled gene expression per time point
  cur.df <- tibble::rownames_to_column(cur.df, "gene")
  h.cond.df <- cur.df %>% left_join(hclusters.df, by = c("gene"))
  colnames(h.cond.df) <- c("gene","0", "6", "12", "24", "48", "53", "cluster")
  h.clustered_expression <- melt(h.cond.df, id = c("gene", "cluster"))
  
  # Generate mean and variance table for each cluster
  h.cluster.means <- h.clustered_expression %>% 
    group_by(cluster, variable) %>% 
    summarize(mean.exp = mean(value, na.rm = TRUE), 
              var = var(value, na.rm = TRUE))
  
  cluster.means.list[[i]] <- h.cluster.means
  
  # Visualize results
  print(h.clustered_expression %>% ggplot(aes(variable, value, group = gene)) 
        + geom_line(alpha = 0.25) 
        + geom_line(aes(variable, 
                        mean.exp, 
                        group = cluster, 
                        color = cluster),
                    data = h.cluster.means,size = 1.1) 
        + ylim(-1.5, 1.5) 
        + facet_wrap(~cluster, ncol = 3) 
        + xlab("Time") 
        + ylab("Scaled Expression") 
        + ggtitle(i))
  
  # Store results
  cluster.info.list[[i]] <- h.cond.df
  
  # Output annotated expression results as csv
  # write.csv(h.cond.df, paste0(i,".gene_clusters.csv"), quote = FALSE)
  
}
