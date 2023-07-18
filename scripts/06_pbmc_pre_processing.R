####################
## Pre-processing ##
####################

# Pre-process and code contributed by Ziyan Lin @Ziyan.Lin@nyulangone.org

##### import necessary libraries #####

#BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Signac)
library(Seurat)
library(ggplot2)
library(GenomicRanges)
library(future)
library(reshape)
library(dplyr)
library(cowplot)
library(stringr)
library(glmGamPoi)

### Input Data ###

sample1 <- Read10X("rawdata/count-0h/")
sample2 <- Read10X("rawdata/count-6h/")
sample3 <- Read10X("rawdata/count-12h/")
sample4 <- Read10X("rawdata/count-24h/")
sample5 <- Read10X("rawdata/count-48h/")
sample6 <- Read10X("rawdata/count-53h/")

### Make Seurat Object ###

sample1.seurat <- CreateSeuratObject(counts = sample1, min.cells = 3, min.features = 200)
sample2.seurat <- CreateSeuratObject(counts = sample2, min.cells = 3, min.features = 200)
sample3.seurat <- CreateSeuratObject(counts = sample3, min.cells = 3, min.features = 200)
sample4.seurat <- CreateSeuratObject(counts = sample4, min.cells = 3, min.features = 200)
sample5.seurat <- CreateSeuratObject(counts = sample5, min.cells = 3, min.features = 200)
sample6.seurat <- CreateSeuratObject(counts = sample6, min.cells = 3, min.features = 200)

# Name the samples

Cond1 <- sample(c("PBMC_2nd_0h"), size = 2500, replace = TRUE)
Cond2 <- sample(c("PBMC_2nd_6h"), size = 8603, replace = TRUE)
Cond3 <- sample(c("PBMC_2nd_12h"), size = 3997, replace = TRUE)
Cond4 <- sample(c("PBMC_2nd_24h"), size = 2539, replace = TRUE)
Cond5 <- sample(c("PBMC_2nd_48h"), size = 8007, replace = TRUE)
Cond6 <- sample(c("PBMC_2nd_53h"), size = 8654, replace = TRUE)

names(Cond1) <- colnames(sample1.seurat)
names(Cond2) <- colnames(sample2.seurat)
names(Cond3) <- colnames(sample3.seurat)
names(Cond4) <- colnames(sample4.seurat)
names(Cond5) <- colnames(sample5.seurat)
names(Cond6) <- colnames(sample6.seurat)

sample1.seurat <- AddMetaData(object = sample1.seurat, metadata = Cond1, col.name = "Condition")
sample2.seurat <- AddMetaData(object = sample2.seurat, metadata = Cond2, col.name = "Condition")
sample3.seurat <- AddMetaData(object = sample3.seurat, metadata = Cond3, col.name = "Condition")
sample4.seurat <- AddMetaData(object = sample4.seurat, metadata = Cond4, col.name = "Condition")
sample5.seurat <- AddMetaData(object = sample5.seurat, metadata = Cond5, col.name = "Condition")
sample6.seurat <- AddMetaData(object = sample6.seurat, metadata = Cond6, col.name = "Condition")


# Merge samples

my.combined <- merge(sample1.seurat, y = c(sample2.seurat,sample3.seurat,sample4.seurat,sample5.seurat,sample6.seurat))


### QC RNA ###

DefaultAssay(my.combined) <- "RNA"
my.combined[["percent.mt"]] <- PercentageFeatureSet(my.combined, pattern = "^MT-")

p <- NULL
p <- VlnPlot(
  object = my.combined,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  group.by = "Condition")
#chartreuse4, chocolate3
png("QC.png",width=10,height=6,res=680,unit="in")
p
dev.off()

#Filtering
my.combined <- subset(my.combined,
                      subset = nFeature_RNA > 200 &
                        percent.mt < 20)
#34300
#31321

#left 91%
p <- NULL
p <- VlnPlot(
  object = my.combined,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  group.by = "Condition")
#chartreuse4, chocolate3
png("QC-norm.png",width=10,height=6,res=680,unit="in")
p
dev.off()


### SCT ###

my.combined <- SCTransform(my.combined, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
DefaultAssay(my.combined) <- "SCT"

# Run the standard workflow for visualization and clustering
my.combined <- ScaleData(my.combined, verbose = FALSE)

# my.combined <- FindVariableFeatures(my.combined)
my.combined <- RunPCA(my.combined, npcs = 30, verbose = FALSE)
my.combined <- RunUMAP(my.combined, reduction = "pca", dims = 1:20)
my.combined <- RunTSNE(my.combined, reduction = "pca", dims = 1:20)
my.combined <- FindNeighbors(my.combined, reduction = "pca", dims = 1:20)
my.combined <- FindClusters(my.combined, resolution = 0.5)

DefaultAssay(my.combined) <- "RNA"
my.combined <- NormalizeData(my.combined)
my.combined <- ScaleData(my.combined, verbose = FALSE)

# Save Seurat Object
saveRDS(my.combined,"my.seurat.rds")

my.combined$Condition <- factor(my.combined$Condition,levels=c("PBMC_2nd_0h","PBMC_2nd_6h","PBMC_2nd_12h","PBMC_2nd_24h","PBMC_2nd_48h","PBMC_2nd_53h"))

# Split Data

split.list <- SplitObject(my.combined, split.by = "Condition")
my.list <- split.list[c("PBMC_2nd_0h","PBMC_2nd_6h","PBMC_2nd_12h","PBMC_2nd_24h","PBMC_2nd_48h","PBMC_2nd_53h")]

### SCT & Integration ###

my.list <- lapply(X = my.list, FUN = SCTransform, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

my.features <- SelectIntegrationFeatures(object.list = my.list, nfeatures = 3000)
my.list <- PrepSCTIntegration(object.list = my.list, anchor.features = my.features)

my.anchors <- FindIntegrationAnchors(object.list = my.list, normalization.method = "SCT",
                                     anchor.features = my.features)

my.combined <- IntegrateData(anchorset = my.anchors, normalization.method = "SCT")

DefaultAssay(my.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
my.combined <- ScaleData(my.combined, verbose = FALSE)
#my.combined <- FindVariableFeatures(my.combined)
my.combined <- RunPCA(my.combined, npcs = 30, verbose = FALSE)
my.combined <- RunUMAP(my.combined, reduction = "pca", dims = 1:20)
my.combined <- RunTSNE(my.combined, reduction = "pca", dims = 1:20)
my.combined <- FindNeighbors(my.combined, reduction = "pca", dims = 1:20)
my.combined <- FindClusters(my.combined, resolution = 0.5)

DefaultAssay(my.combined) <- "RNA"
my.combined <- NormalizeData(my.combined)
my.combined <- ScaleData(my.combined, verbose = FALSE)

my.combined$Condition <- factor(my.combined$Condition,levels=c("PBMC_2nd_0h","PBMC_2nd_6h","PBMC_2nd_12h","PBMC_2nd_24h","PBMC_2nd_48h","PBMC_2nd_53h"))

### Find markers ###

DefaultAssay(my.combined) <- "RNA"

markers <- FindAllMarkers(my.combined)
write.csv(markers,"marker.csv",quote=F)
#markers <- read.csv("marker.csv")

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20,"top20-marker.csv",quote=F)


### Gene plots ###

genes <- readLines("pbmc_genes.txt")

p <- NULL
p <- FeaturePlot(my.combined, features = genes,keep.scale = "all")
png("gene-umap.png",width=16,height=27,res=680,unit="in")
p
dev.off()

#heatmap
# Single cell heatmap of feature expression
heatmap <- DoHeatmap(my.combined , features = genes) + NoLegend()
heatmap2 <- DoHeatmap(subset(my.combined, downsample = 100) , features = genes) + NoLegend()

### Proportion ###

table(my.combined$Condition)
table(my.combined$seurat_clusters)

num <- table(my.combined$seurat_clusters, my.combined$Condition)
write.csv(num,"num.csv",quote=F)
prop <- prop.table(table(my.combined$seurat_clusters, my.combined$Condition), margin = 2)
write.csv(prop,"prop.csv",quote=F)

df <- melt(prop)
colnames(df) <- c("Cluster","Condition","Prop")
df$Condition <- factor(df$Condition,levels=c("Young","ActiveMC"))

p <- NULL
p <- ggplot(df, aes(x=Condition,y=Prop,fill=Condition)) +
  facet_wrap(~Cluster,scales = "free") +
  geom_bar(stat="identity", width=0.7)+
  theme_bw()
p <- p + scale_fill_manual(values=c("tomato","chocolate3"))

#Save Seurat Object
saveRDS(my.combined,"second_PBMC_batchcorr.rds")
