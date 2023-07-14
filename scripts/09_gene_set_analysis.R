###########################################
## Temporally Enriched Gene Set Analysis ##
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
library(reshape2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)

#### Load Clustered Gene Set ####
# If preferred, load output from gene clustering result
# cluster.info.list <- c()
# cluster.means.list <- c()
# 
# for(i in cell_types){
#   cur <- read.csv(paste0(i,".gene_clusters.csv"))
#   cur <- subset(cur, select = -c(X))
#   colnames(cur) <- c("gene","0", "6", "12", "24", "48", "53", "cluster")
#   melted <- melt(cur, id = c("gene", "cluster"))
#   cluster.means <- melted %>% 
#     group_by(variable, cluster) %>% 
#     summarize(mean.exp = mean(value, na.rm = TRUE), 
#               var = var(value, na.rm = TRUE))
#   
#   cluster.info.list[[i]] <- cur
#   cluster.means.list[[i]] <- cluster.means
# }
# 

#### APC MHC II presentation ####

#Note: Of all possible antigen presenting cell types, only the following
#      contained MHC II genes that were temporally variable. Plasma B cells
#      were not included due to low cell count

apc.plots.list <- c()

# List of APC cell type names containing temporally variable MHC II expression
APC <- c( "Monocytes_CD14", 
          "Monocytes_CD16", 
          "Macrophage_1", 
          "Macrophage_2", 
          "Megakaryocyte_Pro")

# List plot names
title_name <- c("mono-CD14", "mono-CD16", "M1", "M2","MGK-P")
id <- 1

for(i in APC){
  
  #load expression info and remove cluster info
  h.cond.df <- cluster.info.list[[i]]
  h.clustered_expression <- melt(h.cond.df, id = c("gene", "cluster"))
  h.clustered_expression["cluster"] <- 
    sample(1, 
           length(h.clustered_expression$cluster), 
           replace = TRUE)
  
  var.genes <- h.cond.df$gene
  
  # Find MHC II genes
  # Note: HLA-D from MHC I is not present in var.genes, 
  #       thus does not require additional removal
  hla <- var.genes[grepl("^HLA-D", var.genes)]
  hla.mean <- h.clustered_expression %>% filter(gene %in% hla) %>% 
    group_by(variable, cluster) %>% 
    summarize(mean.exp = mean(value, na.rm = TRUE), 
              var = var(value, na.rm = TRUE))
  
  # Plot results
  hla.plot <- h.clustered_expression %>% filter(gene %in% hla) %>% 
    ggplot(aes(variable, value, group = gene)) + 
    geom_line(alpha = 0.25) +
    geom_line(aes(variable, mean.exp, group = cluster, color = cluster),
              data = hla.mean,size = 1.1, 
              color = "purple") + 
    theme_classic()+
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          axis.line.x.top = element_blank(),  
          axis.line.y.right = element_blank()) + 
    NoGrid()+ 
    NoLegend() + 
    scale_x_discrete(expand = c(0.05, 0.05)) + 
    ggtitle(title_name[[id]]) + 
    theme(plot.title = element_text(size = 8, hjust = 0.05, vjust = -4.5), 
          axis.text = element_text(size = 8))  +
    ylim(c(-0.5, 0.5))+
    theme(plot.margin = margin(0.3,0.3,0.3,0.3, "mm")) 
  
  id<- id +1
  apc.plots.list[[i]] <- hla.plot
}

## Figure 5.c ##
fig_C <- 
  ggarrange(apc.plots.list$Monocytes_CD14,
            apc.plots.list$Monocytes_CD16,
            apc.plots.list$Macrophage_1,
            apc.plots.list$Macrophage_2,
            apc.plots.list$Megakaryocyte_Pro, nrow = 1, align = "h")

fig_C <- 
  annotate_figure(fig_C, 
                  left = text_grob("Scaled XPR", 
                                   size = 8, 
                                   rot = 90, 
                                   vjust = 0.5), 
                  bottom = text_grob("Hrs Post-Xenotransplantation", 
                                     size = 8))

#################################### Hour 12 ###################################

#### Cluster Selection ####

# Select clusters whose mean expression at 12h > than ave. of mean expression in
# previous time points from each cell type

enriched_clusters12.list <- c()

for(i in cell_types){
  
  # Format loaded gene set info
  diff <- dcast(cluster.means.list[[i]] %>%
                  subset(select = -c(var)), cluster ~ as.character(variable))
  
  # find average of each clusters mean expression from previous time points
  diff$prev <- rowMeans(diff[, c("0", "6")])
  diff$change <- diff[["12"]] - diff[["prev"]]
  
  # filter by difference > 0.5
  diff <- diff %>% filter(change>0.5)
  
  # if cluster exists, record cluster and cell type
  if(length(diff$cluster)>0){
  
    enriched_clusters12.list[[i]] <-diff$cluster
    
  }
}

# Filter selected clusters to only render clusters with >20 genes and >20 cell 
# count for the related cell type for hour 12

# Store plots in list
enrich.12.cluster.plot <- c()

for (i in names(enriched_clusters12.list)){
  
  h.cond.df <- cluster.info.list[[i]]
  h.cluster.expression <- melt(h.cond.df, id = c("gene", "cluster"))
  
  h.cluster.expression <- h.cluster.expression %>% 
    filter(cluster %in% enriched_clusters12.list[[i]])
  
  h.cluster.expression["cluster"] <- 
    sample(1, length(h.cluster.expression$cluster), replace = TRUE)
  
  if(length(h.cluster.expression$gene) > 20){
    
    cell_count <- as.data.frame(table(sample.list[[i]]@meta.data$Time))
    count <- cell_count %>% filter(Var1 == 12)
    count <- count$Freq
    
    if(count>20){
      
      clusters.mean <- h.cluster.expression %>% 
        group_by(variable, cluster) %>% 
        summarize(mean.exp = mean(value, na.rm = TRUE), 
                  var = var(value, na.rm = TRUE))
      
      clusters.mean$cluster <- as.factor(clusters.mean$cluster)
      
      enrich.plot <- h.cluster.expression %>% 
        ggplot(aes(variable, value, group = gene)) + 
        geom_line(alpha = 0.15) +
        geom_line(aes(variable, 
                      mean.exp, 
                      group = cluster, 
                      color = cluster),
                  data = clusters.mean,
                  size = 1.1, color = "skyblue") +
        theme_classic() +
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_blank(), 
              axis.line.x.top = element_blank(),  
              axis.line.y.right = element_blank(),
              axis.text = element_text(size=8)) + 
        NoGrid() + 
        NoLegend() + 
        ylim(c(-1, 1.5)) +
        scale_x_discrete(expand = c(0.05, 0.05)) +
        theme(plot.margin = margin(0.25,0.25,0.25,0.25, "mm"))
      
      enrich.12.cluster.plot[[i]] <- enrich.plot
    }
  }
}

# Cell types as titles for selected clusters belonging to each cell type
e12.names <- c("mCD14","NKT1","M1", "TCD4", "NK","M2", "TCD8", "N. B",
               "TREG", "MGK-P","mCD16","M. B" ,"E-RBC")

# Name plots
id <- 1
for( i in names(enrich.12.cluster.plot)){
  
  enrich.12.cluster.plot[[i]] <- 
    enrich.12.cluster.plot[[i]] + 
    ggtitle(e12.names[[id]]) +
    theme(plot.title = element_text(size = 8, hjust = 0.1, vjust = -4.5))
  
  id <- id +1
}

# Arrange plots
e12.groups <- 
  ggarrange(enrich.12.cluster.plot$Monocytes_CD14,
            enrich.12.cluster.plot$Monocytes_CD16,
            enrich.12.cluster.plot$Macrophage_1,
            enrich.12.cluster.plot$Macrophage_2,
            enrich.12.cluster.plot$NK,
            enrich.12.cluster.plot$NKT1, 
            enrich.12.cluster.plot$T_CD4,
            enrich.12.cluster.plot$T_CD8,
            enrich.12.cluster.plot$T_Reg,
            enrich.12.cluster.plot$Naive_B,
            enrich.12.cluster.plot$Memory_B,
            enrich.12.cluster.plot$Megakaryocyte_Pro,
            ncol = 3, nrow = 4, align = "hv") +
  theme(plot.margin = margin(0,0,0,0, "mm"))

#### 12H Gene Set Expression ####

## Figure 5.d ##
fig_D <- 
  annotate_figure(e12.groups,
                  left = text_grob("Scaled Expression", rot = 90, size = 9),
                  bottom = text_grob("Hours Post-Xenotransplantation",size = 9))

#### 12H Selected Gene Sets Average Normalized Expression ####

# For each cell type included in selected clusters:
# Rerun preprocessing (normalize by feature on corrected UMI, scaling)
# Store in temporary Seurat object, add time point meta data
ave12.df <- data.frame()

for(i in names(enrich.12.cluster.plot)){
  
  cur_cell_type <- sample.list[[i]]
  cur_genes_20 <- genes.list[[i]]
  mtx <- GetAssayData(cur_cell_type, assay = "SCT", slot = "count")
  data <- CreateSeuratObject(counts = mtx, assay = "SCT")
  data <- AddMetaData(data, cur_cell_type@meta.data$Time, col.name = "Time")
  DefaultAssay(data) <- "SCT"
  data <- NormalizeData(data, normalization.method = "CLR", margin = 1)
  data <- ScaleData(data, assay = "SCT", features = rownames(data))
  Idents(data) <- data@meta.data$Time
  data@meta.data$Time <- factor( data@meta.data$Time,
                                 levels = c("0", "6", "12", "24", "48", "53"))
  
  # Extract genes from selected gene set
  genes <- as.data.frame(enrich.12.cluster.plot[[i]]$data)
  genes <- genes$gene
  
  # Average normalized expression values by time point
  ave <- AverageExpression(data, assay = "SCT", 
                           slot = "data", 
                           features = genes, 
                           group.by = "Time")
  
  ave <- as.data.frame(ave$SCT)
  
  # Annotate cell type
  ave["cell_type"] <- sample(i, length(rownames(ave)), replace = TRUE)
  
  # Combine results
  ave12.df <- rbind(ave12.df, ave)
  
}

## Supplementary Table E12 ##
ave12 <- tibble::rownames_to_column(ave12.df, "gene")

write.csv(ave12,"Xeno-Supplementary-Table-E12.csv", quote = FALSE)

########## Hour 12 GO ##########

#### BULK ####

# Combine genes from selected gene sets
combined_genes12 <- c()

for (i in names(enrich.12.cluster.plot)){
  
  peak.clusters <- enriched_clusters12.list[[i]]
  cluster.info <- cluster.info.list[[i]]
  peak.cluster.info <- cluster.info %>% filter(cluster %in% peak.clusters)
  genes <- peak.cluster.info$gene
  combined_genes12 <- c(combined_genes12, genes)
  
}

# Drop duplicates
combined12 <- unique(combined_genes12)

# Extract backgrounds genes
Idents(sample) <- sample@meta.data$cell_type
sub12 <- subset(sample, ident = c(names(enrich.12.cluster.plot)))

DefaultAssay(sample) <- "RNA"
bg12 <- AverageExpression(sub12, 
                          assay = "RNA", 
                          slot = "data", 
                          features = rownames(sub12), 
                          group.by = "orig.ident")
bg12 <- as.data.frame(bg12$RNA)
bg12.filt <- bg12 %>% filter(all>0)

# Perform GO term analysis
enrich.res12 <- enrichGO(gene = combined12,
                         OrgDb = org.Hs.eg.db, 
                         keyType = 'SYMBOL',
                         ont = "ALL", 
                         pAdjustMethod = "bonferroni",
                         universe = rownames(bg12.filt),
                         pvalueCutoff = 0.01, 
                         qvalueCutoff = 0.99, 
                         pool = TRUE)

# Filter GO terms by p.adjust values and gene count
enrich.res12.df <- enrich.res12@result %>% 
  filter(p.adjust <= 5e-07) %>% 
  filter(Count > 15)

# Extract top 5 terms whose gene sets are variable
e12.top5 <- enrich.res12.df %>% 
  filter(ID == "GO:0042110"|
           ID == "GO:0007249" |
           ID == "GO:0097193" |
           ID == "GO:0009896" | 
           ID == "GO:1903131")

#### 12H Bulk GO Term Plot ####

## Figure 5.f ##
fig_F <-
  e12.top5 %>% ggplot(aes(x = reorder(Description, -log10(p.adjust)), 
                          y =-log10(p.adjust), 
                          fill = ONTOLOGY)) +
  geom_bar(stat = "identity", 
           width = 0.5, 
           position = "dodge", 
           fill = "skyblue") +
  theme_classic() + 
  ggtitle("Hour 12 Enrichment") + 
  ylab("-log10(p.adjust)") +
  theme(axis.text = element_text(size = 9), 
        axis.title = element_text(size = 9),
        plot.title = element_text(size = 9, vjust = -1.25),
        axis.title.y = element_blank())+ 
  coord_flip() + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 37))+
  NoLegend()

fig_F <- fig_F + ggtitle("12H Enrichment") + 
  theme(plot.title = element_text(size = 9)) +
  scale_y_continuous(breaks = c(0,3,6,9,12), limits = c(0, 12)) 

#### CELL TYPE SPECIFIC ####

# Extract selected clusters for each cell type
# Use all genes with detected expression for that cell type as background genes
# Perform GO term analysis
cell_type12.res <- c()

for(i in names(enrich.12.cluster.plot)){
  
  input_genes <- unique(enrich.12.cluster.plot[[i]]$data$gene)
  cur <- sample.list[[i]]
  Idents(cur) <- cur@meta.data$cell_type
  
  backg.genes <- AverageExpression(sample.list[[i]], 
                                   assay = "RNA", 
                                   slot = "data")
  
  backg.genes <- as.data.frame(backg.genes$RNA)
  backg.genes <- backg.genes %>% filter(all > 0)
  res <- enrichGO(gene = input_genes,
                  OrgDb = org.Hs.eg.db, 
                  keyType = 'SYMBOL',
                  ont = "ALL", 
                  pAdjustMethod = "bonferroni",
                  universe = rownames(backg.genes),
                  pvalueCutoff = 0.01, 
                  qvalueCutoff = 0.99, 
                  pool = TRUE)
  
  cell_type12.res[[i]] <- res
  
}

# Filter GO terms by p.adjust and gene count and combine results
selected.go12 <- data.frame()

for(i in names(cell_type12.res)){
  
  cur.df <- cell_type12.res[[i]]@result
  cur.df <- cur.df %>% filter(p.adjust < 1e-3)
  cur.df <- cur.df %>% filter(Count > 5)
  cur.df["cell_type"] <- sample(i, length(cur.df$ID), replace = TRUE)
  selected.go12 <- rbind(selected.go12, cur.df)
  
}

# Select top 5 go terms with variable gene sets, at least 1 per cell type
cell_type.e12.top5 <- selected.go12 %>% 
  filter(ID == "GO:0005802" |
           ID == "GO:0005667" |
           ID == "GO:0050727" |
           ID == "GO:2001234" |
           ID == "GO:0031331")

# Annotate by cell type
cell_type.e12.top5["cell_type"] <- 
  c("mono-CD14", "NKT1", "NK", "M2", "mono-CD16")

#### 12H Cell Type Specific GO Term Plot ####

## Supplementary Figure 7 Pt.1 ##
sup.12.go <- cell_type.e12.top5 %>% 
  ggplot(aes(x = reorder(Description, -log10(p.adjust)), 
             y =-log10(p.adjust), 
             fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.5, position = "dodge") +
  theme_classic() + 
  ggtitle("Hour 12 Enrichment") + 
  ylab("-log10(p.adjust)")+
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 8, vjust = -1.25),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_blank())+ 
  coord_flip() + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))

#################################### Hour 48 ###################################

#### Cluster Selection ####

# Select clusters whose mean expression at 48h > than ave. of mean expression in
# previous time points from each cell type

enriched_clusters48.list <- c()

for(i in cell_types){
  
  diff <- dcast(cluster.means.list[[i]] %>%
                  subset(select = -c(var)), cluster ~ as.character(variable))
  
  diff$prev <- rowMeans(diff[, c("0", "6", "12", "24")])
  diff$change <- diff[["48"]] - diff[["prev"]]
  diff <- diff %>% filter(change> 0.50)
  
  if(length(diff$cluster)>0){
    
    enriched_clusters48.list[[i]] <-diff$cluster
    
  }
}

# Filter selected clusters to only render clusters with >20 genes and >20 cell 
# count for the related cell type for hour 12
enrich.48.cluster.plot <- c()
e48.df <- c()

for (i in names(enriched_clusters48.list)){
  
  h.cond.df <- cluster.info.list[[i]]
  h.cluster.expression <- melt(h.cond.df, id = c("gene", "cluster"))
  h.cluster.expression <- h.cluster.expression %>% 
    filter(cluster %in% enriched_clusters48.list[[i]])
  
  h.cluster.expression["cluster"] <- 
    sample(1, length(h.cluster.expression$cluster), replace = TRUE)
  
  if(length(h.cluster.expression$gene) > 20){
    
    cell_count <- as.data.frame(table(sample.list[[i]]@meta.data$Time))
    count <- cell_count %>% filter(Var1 == 48)
    count <- count$Freq
    
    if(count>20){
      
      clusters.mean <- h.cluster.expression %>% 
        group_by(variable, cluster) %>% 
        summarize(mean.exp = mean(value, na.rm = TRUE), 
                  var = var(value, na.rm = TRUE))
      
      clusters.mean$cluster <- as.factor(clusters.mean$cluster)
      
      enrich.plot <- h.cluster.expression %>% 
        ggplot(aes(variable, value, group = gene)) + 
        geom_line(alpha = 0.15) +
        geom_line(aes(variable, 
                      mean.exp, 
                      group = cluster, 
                      color = cluster), 
                  data = clusters.mean,size = 1.1, 
                  color = "tomato") + 
        theme_classic()+
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_blank(), 
              axis.line.x.top = element_blank(),  
              axis.line.y.right = element_blank(),
              axis.text = element_text(size=8)) + 
        NoGrid() + 
        NoLegend() + 
        ylim(c(-1,1.5)) +
        scale_x_discrete(expand = c(0.05, 0.05)) +
        theme(plot.margin = margin(0.25,0.25,0.25,0.25, "mm"))
      
      enrich.48.cluster.plot[[i]] <- enrich.plot
      e48.df[[i]] <- c(h.cluster.expression$gene)
      
    }
  }
}

# Cell types as titles for selected clusters belonging to each cell type
e48.names <- c("NKT1", "MGK","TCD4","NK","M2", "TReg","MGK-P","E-RBC")

# Name plots
id <- 1
for( i in names(enrich.48.cluster.plot)){
  enrich.48.cluster.plot[[i]] <- 
    enrich.48.cluster.plot[[i]] + 
    ggtitle(e48.names[[id]]) +
    theme(plot.title = element_text(size = 8, hjust = 0.1, vjust = -4.5))
  
  id <- id +1
}

# Arrange plots
e48.groups <- 
  ggarrange(enrich.48.cluster.plot$Macrophage_2,
            enrich.48.cluster.plot$NK,
            enrich.48.cluster.plot$NKT1,
            enrich.48.cluster.plot$T_CD4,
            enrich.48.cluster.plot$T_Reg,
            enrich.48.cluster.plot$Megakaryocyte,
            enrich.48.cluster.plot$Megakaryocyte_Pro,
            ncol = 2, nrow = 4, align = "hv") +
  theme(plot.margin = margin(0,0,0,0, "mm"))

#### 48H Gene Set Expression ####

## Figure5.e ##
fig_E<- 
  annotate_figure(e48.groups, left = text_grob("Scaled Expression", 
                                               rot = 90, 
                                               size = 8, 
                                               hjust = -0.55), 
                  bottom = text_grob("Hours Post-Xenotransplantation", 
                                     size = 8, 
                                     vjust = -11))

#### 48H Selected Gene Sets Average Normalized Expression ####
ave48.df <- data.frame()

# For each cell type included in selected clusters:
# Rerun preprocessing (normalize by feature on corrected UMI, scaling)
# Store in temporary Seurat object, add time point meta data
for(i in names(enrich.48.cluster.plot)){
  
  cur_cell_type <- sample.list[[i]]
  cur_genes_20 <- genes.list[[i]]
  mtx <- GetAssayData(cur_cell_type, assay = "SCT", slot = "count")
  data <- CreateSeuratObject(counts = mtx, assay = "SCT")
  data <- AddMetaData(data, cur_cell_type@meta.data$Time, col.name = "Time")
  DefaultAssay(data) <- "SCT"
  data <- NormalizeData(data, normalization.method = "CLR", margin = 1)
  data <- ScaleData(data, assay = "SCT", features = rownames(data))
  Idents(data) <- data@meta.data$Time
  data@meta.data$Time <- factor( data@meta.data$Time,
                                 levels = c("0", "6", "12", "24", "48", "53"))
  
  # Extract genes within selected gene set
  genes <- as.data.frame(enrich.48.cluster.plot[[i]]$data)
  genes <- genes$gene
  
  # Average normalized expression values by time point
  ave <- AverageExpression(data, 
                           assay = "SCT", 
                           slot = "data", 
                           features = genes, 
                           group.by = "Time")
  
  ave <- as.data.frame(ave$SCT)
  
  # Annotate cell type
  ave["cell_type"] <- sample(i, length(rownames(ave)), replace = TRUE)
  
  # Combine results
  ave48.df <- rbind(ave48.df, ave)
  
}

## Supplementary Table E48 ##
ave48 <- tibble::rownames_to_column(ave48.df, "gene")

write.csv(ave48,"Xeno-Supplementary-Table-E48.csv", quote = FALSE)

########## Hour 48 GO ##########

#### BULK #### 

# Combine genes from selected gene sets
combined_genes48 <- c()

for (i in names(enrich.48.cluster.plot)){
  peak.clusters <- enriched_clusters48.list[[i]]
  
  cluster.info <- cluster.info.list[[i]]
  
  peak.cluster.info <- cluster.info %>% filter(cluster %in% peak.clusters)
  
  genes <- peak.cluster.info$gene
  
  combined_genes48 <- c(combined_genes48, genes)
}

# Drop duplicates
combined48 <- unique(combined_genes48)

# Extract background genes
Idents(sample) <- sample@meta.data$cell_type
sub48 <- subset(sample, ident = c(names(enrich.48.cluster.plot)))

DefaultAssay(sample) <- "RNA"
bg48 <- AverageExpression(sub48, 
                          assay = "RNA", 
                          slot = "data", 
                          features = rownames(sub48), 
                          group.by = "orig.ident")

bg48 <- as.data.frame(bg48$RNA)
bg48.filt <- bg48 %>% filter(all>0)

# Perform GO term analysis
enrich.res48 <- enrichGO(gene = combined48,
                         OrgDb = org.Hs.eg.db, 
                         keyType = 'SYMBOL',
                         ont = "ALL", 
                         pAdjustMethod = "bonferroni",
                         universe = rownames(bg48.filt),
                         pvalueCutoff = 0.01, 
                         qvalueCutoff = 0.99, 
                         pool = TRUE)

# Filter GO terms by p.adjust values and gene count
enrich.res48.df <- enrich.res48@result %>% 
  filter(p.adjust <= 5e-07)%>% 
  filter(Count > 5)

# Extract top 5 terms whose gene sets are variable
e48.top5 <- enrich.res48.df %>% filter(ID == "GO:0042110"|
                                             ID == "GO:0001906" | 
                                             ID == "GO:0034341" |
                                             ID == "GO:0032943" | 
                                             ID == "GO:0032103")

#### 48H Bulk GO Term Plot ####

## Figure 5.g ##
fig_G <- 
  e48.top5 %>% ggplot(aes(x = reorder(Description, -log10(p.adjust)), 
                                y =-log10(p.adjust), 
                                fill = ONTOLOGY)) +
  geom_bar(stat = "identity", 
           width = 0.5, 
           position = "dodge", 
           fill = "tomato") +
  theme_classic() + 
  ggtitle("Hour 48 Enrichment") + 
  ylab("-log10(p.adjust)") +
  theme(axis.text = element_text(size = 9), 
        axis.title = element_text(size = 9),
        plot.title = element_text(size = 8, vjust = -1.25),
        axis.title.y = element_blank()) + 
  coord_flip() + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
  NoLegend()

fig_G <- fig_G +
  ggtitle("48H Enrichment") + 
  theme(plot.title = element_text(size = 9))

#### CELL TYPE SPECIFIC #### 

# Extract selected clusters for each cell type
# Use all genes with detected expression for that cell type as background genes
# Perform GO term analysis
cell_type48.res <- c()

for(i in names(enrich.48.cluster.plot)){
  
  input_genes <- unique(enrich.48.cluster.plot[[i]]$data$gene)
  
  cur <- sample.list[[i]]
  
  Idents(cur) <- cur@meta.data$cell_type
  
  backg.genes <- AverageExpression(sample.list[[i]], assay = "RNA", slot = "data")
  
  backg.genes <- as.data.frame(backg.genes$RNA)
  
  backg.genes <- backg.genes %>% filter(all > 0)
  
  res <- enrichGO(gene = input_genes,
                  OrgDb = org.Hs.eg.db, 
                  keyType = 'SYMBOL',
                  ont = "ALL", 
                  pAdjustMethod = "bonferroni",
                  universe = rownames(backg.genes),
                  pvalueCutoff = 0.01, 
                  qvalueCutoff = 0.99, 
                  pool = TRUE)
  
  cell_type48.res[[i]] <- res
  
}

# Filter GO terms by p.adjust and gene count and combine results
selected.go48 <- data.frame()

for(i in names(cell_type48.res)){
  
  cur.df <- cell_type48.res[[i]]@result
  
  cur.df <- cur.df %>% filter(p.adjust < 1e-3)
  
  cur.df <- cur.df %>% filter(Count > 5)
  
  cur.df["cell_type"] <- sample(i, length(cur.df$ID), replace = TRUE)
  
  selected.go48 <- rbind(selected.go48, cur.df)
  
}

# Select top 5 go terms, at least 1 by cell type, with variable gene sets
e48.nkt1 <- selected.go48 %>% filter(cell_type == "NKT1") %>% 
  filter(ID == "GO:0044194")

e48.tcd4 <- selected.go48 %>% filter(cell_type == "T_CD4") %>% 
  filter(ID == "GO:0034341" |ID == "GO:0042110"|
           ID == "GO:0051249"| ID == "GO:0050852") 

e48.nk <- selected.go48 %>% filter(cell_type == "NK") %>%
  filter(ID == "GO:0034341")

e48.treg <- selected.go48 %>% filter(cell_type == "T_Reg") %>%
  filter(ID == "GO:0034341" |ID == "GO:0042110" | ID == "GO:0050852")

# Combine
cell_type.e48.top5 <- rbind(e48.nkt1, e48.nk, e48.tcd4, e48.treg)

# Annotate by cell type
cell_type.e48.top5["cell_type"] <- c("NKT1", "NK", "T-CD4", "T-CD4","T-CD4","T-CD4","T-Reg","T-Reg","T-Reg")

#### 48H Cell Type Specific GO Term Plot ####

## Supplementary Figure 7 Pt.2 ##
sup.48.go <- cell_type.e48.top5 %>% 
  ggplot(aes(x = reorder(Description, -log10(p.adjust)), 
             y =-log10(p.adjust), 
             fill = cell_type)) +
  geom_bar(stat = "identity", 
           position = position_dodge2(preserve = "single")) +
  theme_classic() + 
  ggtitle("Hour 48 Enrichment") + 
  ylab("-log10(p.adjust)")+
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 8, vjust = -1.25),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_blank())+ 
  coord_flip() + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))

#### Cell Type Specific GO Term Plot ####

## Supplementary Figure 7 ##
fig_S7 <- ggarrange(sup.12.go, sup.48.go, nrow = 1, align = "h")