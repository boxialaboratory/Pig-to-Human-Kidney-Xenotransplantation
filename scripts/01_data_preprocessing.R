##########################
##### Pre-processing #####
##########################

##### import libraries #####
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(rhdf5)

##### Step 1: Find Homologous Genes between Human and Pig #####
# Starting database is downloadable from: https://download.cncb.ac.cn/hgd/homolog_gene/Pig_Human_Homolog_gene.txt.gz
# This file is not included in the data folder as its size 121.21MB exceeds Github's maximum capacity per file (100MB).
# An edited version is included in the data folder. Code that performed this edit is below but commented out.

# homolog_database <- read.table(file = "data/Pig_Human_Homolog_gene_OG.txt", header = T, fill = T, row.names = NULL)
# homolog_database <- homolog_database[, c(4, 6, 12, 14)] # 554797 genes

# # rename the columns
# colnames(homolog_database) <- c("organism1", "gene1", "organism2", "gene2")

# # Filter the rows containing "Human" or "Pig" in organism1 and organism2. Other values include eggNOG:xxx; TreeFam: xxx; Panther:xxx
# homolog_filtered <- homolog_database %>%
#   filter((organism1 == "Human" & organism2 == "Pig") | (organism2 == "Human" & organism1 == "Pig"))

# # 529988 homo-genes left. Store edited version of database.
# write.csv(homolog_filtered, "data/homolog_filtered.csv")

# Load in edited version of database
homolog_filtered <- read.csv("data/homolog_filtered.csv")

homolog_filtered <- homolog_filtered %>%
  filter(gene1 == gene2) %>%
  distinct(gene1, gene2, .keep_all = TRUE)
# 14519 genes
homolog_list = homolog_filtered$gene1 # human_genes = pig_genes

##### Step 2: Create Expression Matrices #####
# Kidney 1:
# note: A and B are technical replicates of each other
mat.human.t.A <- Read10X_h5("data/first_kidney_transplanted_A_hg38_filtered_feature_bc_matrix.h5")
mat.human.t.B <- Read10X_h5("data/first_kidney_transplanted_B_hg38_filtered_feature_bc_matrix.h5")
mat.human.u.A <- Read10X_h5("data/first_kidney_untransplanted_A_hg38_filtered_feature_bc_matrix.h5")
mat.human.u.B <- Read10X_h5("data/first_kidney_untransplanted_B_hg38_filtered_feature_bc_matrix.h5")
# add "_1" and "_2" to cell names for human untreated condition to align with other conditions
mat.human.u.A@Dimnames[[2]] <- paste(mat.human.u.A@Dimnames[[2]],"1",sep="_")
mat.human.u.B@Dimnames[[2]] <- paste(mat.human.u.B@Dimnames[[2]],"2",sep="_")

# mapped to porcine
mat.pig.t.A <- Read10X_h5("data/first_kidney_transplanted_A_Sscrofa11_filtered_feature_bc_matrix.h5")
mat.pig.t.B <- Read10X_h5("data/first_kidney_transplanted_B_Sscrofa11_filtered_feature_bc_matrix.h5")
mat.pig.u.A <- Read10X_h5("data/first_kidney_untransplanted_A_Sscrofa11_filtered_feature_bc_matrix.h5")
mat.pig.u.B <- Read10X_h5("data/first_kidney_untransplanted_B_Sscrofa11_filtered_feature_bc_matrix.h5")

#make seurat object
human.t.A <- CreateSeuratObject(counts = mat.human.t.A)
human.t.B <- CreateSeuratObject(counts = mat.human.t.B)
human.u.A <- CreateSeuratObject(counts = mat.human.u.A)
human.u.B <- CreateSeuratObject(counts = mat.human.u.B)

pig.t.A <- CreateSeuratObject(counts = mat.pig.t.A)
pig.t.B <- CreateSeuratObject(counts = mat.pig.t.B)
pig.u.A <- CreateSeuratObject(counts = mat.pig.u.A)
pig.u.B <- CreateSeuratObject(counts = mat.pig.u.B)

# merge technical duplicates into one object
human.t <- merge(human.t.A, human.t.B)
human.u <- merge(human.u.A, human.u.B)
pig.t <- merge(pig.t.A, pig.t.B)
pig.u <- merge(pig.u.A, pig.u.B)

# there are some duplicates cells so the barcodes are renamed to ensure uniqueness
# Add Conditions as metadata
Cond1 <- sample(c("kidney_1st_treated"), size = 6301, replace = TRUE)
Cond2 <- sample(c("kidney_1st_untreated_rt"), size = 2436, replace = TRUE)

names(Cond1) <- colnames(human.t)
names(Cond2) <- colnames(human.u)

human.t <- AddMetaData(object = human.t, metadata = Cond1, col.name = "Condition")
human.u <- AddMetaData(object = human.u, metadata = Cond2, col.name = "Condition")

Cond1 <- sample(c("kidney_1st_treated"), size = 9019, replace = TRUE)
Cond2 <- sample(c("kidney_1st_untreated_rt"), size = 4616, replace = TRUE)

names(Cond1) <- colnames(pig.t)
names(Cond2) <- colnames(pig.u)

pig.t <- AddMetaData(object = pig.t, metadata = Cond1, col.name = "Condition")
pig.u <- AddMetaData(object = pig.u, metadata = Cond2, col.name = "Condition")

##### Step 3: Quality Control #####
# QC: genes detected > 500; mitochondrial genes < 8%; UMI > 1000
qc <- function(data) {
  DefaultAssay(data) <- "RNA"
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  #Filtering
  data <- subset(data, subset = nFeature_RNA > 500 & percent.mt < 8 & nCount_RNA > 1000)
  return(data)
}

# note that mitochondrial genes were not annotated in the porcine-mapped dataset
# all pig-mapped cells will show having 0 "percent.mt"

human.t <- qc(human.t)
human.u <- qc(human.u)
pig.t <- qc(pig.t)
pig.u <- qc(pig.u)

# visualize distribution of quality control metrics
VlnPlot(human.t, # human.u, pig.t, pig.u
        features = c("nFeature_RNA","nCount_RNA","percent.mt"),
        ncol = 3,pt.size = 0.1) & theme(plot.title = element_text(size=10))

##### Step 4: Create Homo-gene only Human and Pig Matrices #####
human.t.homolog <- subset(human.t,
                          features = homolog_list)
human.u.homolog <- subset(human.u,
                          features = homolog_list)
pig.t.homolog <- subset(pig.t,
                        features = homolog_list)
pig.u.homolog <- subset(pig.u,
                        features = homolog_list)

### Step 5: UMI MAPPING ###
# Treated condition #
# find all cells that are in human and not in pig
unique_human_cells <- setdiff(colnames(human.t.homolog@assays$RNA@counts), colnames(pig.t.homolog@assays$RNA@counts)) 

# find all cells that are in pig and not in human
unique_porcine_cells <- setdiff(colnames(pig.t.homolog@assays$RNA@counts), colnames(human.t.homolog@assays$RNA@counts))

# find all the cells that are duplicated in both
duplicated_cells <- intersect(colnames(pig.t.homolog@assays$RNA@counts), colnames(human.t.homolog@assays$RNA@counts))

# Subset the human and pig Seurat objects to only contain the unique cells
human.t.homolog.unique <- human.t.homolog[, unique_human_cells]
pig.t.homolog.unique <- pig.t.homolog[, unique_porcine_cells]

# Merge the Seurat objects so duplicated cells are only included once (union)
merged.t <- merge(human.t.homolog.unique, 
                  y = pig.t.homolog, 
                  project = "Merged")

# Add "origin" metadata to the seurat objects with their host organism assignment
# for human-only and porcine-only cells, assign the metadata accordingly
merged.t@meta.data$origin <- NULL
merged.t@meta.data$percent_umi_human <- NULL
merged.t@meta.data$percent_umi_porcine <- NULL
merged.t@meta.data$origin[colnames(merged.t@assays$RNA@counts) %in% unique_human_cells] <- "human"
merged.t@meta.data$origin[colnames(merged.t@assays$RNA@counts) %in% unique_porcine_cells] <- "porcine"

# For unique human cells, set percent_umi_human to 100 and percent_umi_porcine to 0
merged.t@meta.data$percent_umi_human[colnames(merged.t) %in% unique_human_cells] <- 100
merged.t@meta.data$percent_umi_porcine[colnames(merged.t) %in% unique_human_cells] <- 0
merged.t@meta.data$percent_umi_porcine[colnames(merged.t) %in% unique_porcine_cells] <- 100
merged.t@meta.data$percent_umi_human[colnames(merged.t) %in% unique_porcine_cells] <- 0

human_umi_counts <- Matrix::colSums(human.t.homolog@assays$RNA@counts)
porcine_umi_counts <- Matrix::colSums(pig.t.homolog@assays$RNA@counts)

for (cell in duplicated_cells) {
  # fill in "origin" metadata with assigned organism
  if (human_umi_counts[cell] > porcine_umi_counts[cell]) {
    merged.t@meta.data$origin[colnames(merged.t@assays$RNA@counts) == cell] <- "human"
  } else {
    merged.t@meta.data$origin[colnames(merged.t@assays$RNA@counts) == cell] <- "porcine"
    # fill in percent_umi_human with percentage of umi mapped to human
  }
    total_umi <- human_umi_counts[cell] + porcine_umi_counts[cell]
    merged.t@meta.data$percent_umi_human[colnames(merged.t) == cell] <- round((human_umi_counts[cell] / total_umi[cell]) * 100, 2)
    # fill in percent_umi_human with percentage of umi mapped to porcine
    merged.t@meta.data$percent_umi_porcine[colnames(merged.t) == cell] <- round((porcine_umi_counts[cell] / total_umi[cell]) * 100, 2)
  }

table(merged.t@meta.data$origin)

# Look at UMI density
# convert the data to a data frame
df <- data.frame(
  percent_umi_human = merged.t@meta.data$percent_umi_human,
  percent_umi_porcine = merged.t@meta.data$percent_umi_porcine
)

# Add a new column 'color' to the dataframe based on the condition
df$color <- ifelse(df$percent_umi_porcine > 50, "red", "blue")

# Generate the plot
ggplot(df, aes(x = percent_umi_human, y = percent_umi_porcine, color = color)) +
  geom_point(alpha = 0.5) +
  scale_color_identity() +
  theme_minimal() +
  labs(
    x = "Percent UMI Human",
    y = "Percent UMI Porcine",
    title = "Treated Kidney 1 Percent UMI Human vs Percent UMI Porcine") + theme_bw() + ylim(0,100) + xlim(0,100)

### create matrix with human and pig cells with assigned origin ###
human_cell_barcodes <- merged.t@assays$RNA@counts@Dimnames[[2]][merged.t@meta.data$origin == "human"]
human.cells <- human.t.homolog[, human_cell_barcodes]
human.cells@meta.data$origin <- rep("human", nrow(human.cells@meta.data))

# Subset pig.t to exclude cells with matching barcodes
pig_cell_barcodes <- merged.t@assays$RNA@counts@Dimnames[[2]][merged.t@meta.data$origin == "porcine"]
pig.cells <- pig.t.homolog[, pig_cell_barcodes]
pig.cells@meta.data$origin <- rep("porcine", nrow(pig.cells@meta.data))

treated <- merge(human.cells, 
                  pig.cells, 
                  project = "Merged")

# Untreated #
# Find the unique and duplicated cells between the datasets
# find all cells that are in human and not in pig
unique_human_cells <- setdiff(colnames(human.u.homolog@assays$RNA@counts), colnames(pig.u.homolog@assays$RNA@counts)) # 0 unique human cells

# find all cells that are in pig and not in human
unique_porcine_cells <- setdiff(colnames(pig.u.homolog@assays$RNA@counts), colnames(human.u.homolog@assays$RNA@counts)) # 1144 cells

# find all the cells that are duplicated in both
duplicated_cells <- intersect(colnames(pig.u.homolog@assays$RNA@counts), colnames(human.u.homolog@assays$RNA@counts)) # 1675 cells

# All the cell barcodes for human-mapped are contained within the porcine matrix
merged.u <- pig.u.homolog

# Add "origin" metadata to the seurat objects with their host organism assignment
# for human-only and porcine-only cells, assign the metadata accordingly
merged.u@meta.data$origin <- NULL
merged.u@meta.data$percent_umi_human <- NULL
merged.u@meta.data$percent_umi_porcine <- NULL
merged.u@meta.data$origin[colnames(merged.u@assays$RNA@counts) %in% unique_porcine_cells] <- "porcine"

# For unique human cells, set percent_umi_human to 100 and percent_umi_porcine to 0
merged.u@meta.data$percent_umi_human[colnames(merged.u) %in% unique_human_cells] <- 100
merged.u@meta.data$percent_umi_porcine[colnames(merged.u) %in% unique_human_cells] <- 0
merged.u@meta.data$percent_umi_porcine[colnames(merged.u) %in% unique_porcine_cells] <- 100
merged.u@meta.data$percent_umi_human[colnames(merged.u) %in% unique_porcine_cells] <- 0

human_umi_counts <- Matrix::colSums(human.u.homolog@assays$RNA@counts)
porcine_umi_counts <- Matrix::colSums(pig.u.homolog@assays$RNA@counts)

for (cell in duplicated_cells) {
  # fill in "origin" metadata with assigned organism
  if (human_umi_counts[cell] > porcine_umi_counts[cell]) {
    merged.u@meta.data$origin[colnames(merged.u@assays$RNA@counts) == cell] <- "human"
  } else {
    merged.u@meta.data$origin[colnames(merged.u@assays$RNA@counts) == cell] <- "porcine"
  }
    # fill in percent_umi_human with percentage of umi mapped to human
    total_umi <- human_umi_counts[cell] + porcine_umi_counts[cell]
    merged.u@meta.data$percent_umi_human[colnames(merged.u@assays$RNA@counts) == cell] <- round((human_umi_counts[cell] / total_umi) * 100, 2)
    # fill in percent_umi_human with percentage of umi mapped to porcine
    merged.u@meta.data$percent_umi_porcine[colnames(merged.u@assays$RNA@counts) == cell] <- round((porcine_umi_counts[cell] / total_umi) * 100, 2)
  }

table(merged.u@meta.data$origin)
# porcine 
# 2819

# Look at UMI density
# convert the data to a data frame
df <- data.frame(
  percent_umi_human = merged.u@meta.data$percent_umi_human,
  percent_umi_porcine = merged.u@meta.data$percent_umi_porcine
)
# Add a new column 'color' to the dataframe based on the condition
df$color <- ifelse(df$percent_umi_porcine > 50, "red", "blue")

# Generate the plot
ggplot(df, aes(x = percent_umi_human, y = percent_umi_porcine, color = color)) +
  geom_point(alpha = 0.5) +
  scale_color_identity() +
  theme_minimal() +
  labs(
    x = "Percent UMI Human",
    y = "Percent UMI Porcine",
    title = "Untreated Kidney 1 Percent UMI Human vs Percent UMI Porcine") + theme_bw() + ylim(0,100) + xlim(0,100)

untreated <- merged.u # rename for organizational purposes

#########################################################################
######################### repeat with kidney 2 ########################## 
#########################################################################

##### Step 2: Create Expression Matrices #####
mat.human.t <- Read10X_h5("data/second_kidney_transplanted_hg38_filtered_feature_bc_matrix.h5")
mat.human.urt <- Read10X_h5("data/second_kidney_untransplanted_rt_hg38_filtered_feature_bc_matrix.h5")

# mapped to porcine
mat.pig.t <- Read10X_h5("data/second_kidney_transplanted_Sscrofa11_filtered_feature_bc_matrix.h5")
mat.pig.urt <- Read10X_h5("data/second_kidney_untransplanted_rt_Ssacrofa11_filtered_feature_bc_matrix.h5")

#make seurat object
human.t <- CreateSeuratObject(counts = mat.human.t)
human.urt <- CreateSeuratObject(counts = mat.human.urt)

pig.t <- CreateSeuratObject(counts = mat.pig.t)
pig.urt <- CreateSeuratObject(counts = mat.pig.urt)

# Name the samples
Cond1 <- sample(c("kidney_2nd_treated"), size = 4814, replace = TRUE)
Cond2 <- sample(c("kidney_2nd_untreated_rt"), size = 7032, replace = TRUE)

names(Cond1) <- colnames(human.t)
names(Cond2) <- colnames(human.urt)

human.t <- AddMetaData(object = human.t, metadata = Cond1, col.name = "Condition")
human.urt <- AddMetaData(object = human.urt, metadata = Cond2, col.name = "Condition")

Cond1 <- sample(c("kidney_2nd_treated"), size = 7291, replace = TRUE)
Cond2 <- sample(c("kidney_2nd_untreated_rt"), size = 10755, replace = TRUE)

names(Cond1) <- colnames(pig.t)
names(Cond2) <- colnames(pig.urt)

pig.t <- AddMetaData(object = pig.t, metadata = Cond1, col.name = "Condition")
pig.urt <- AddMetaData(object = pig.urt, metadata = Cond2, col.name = "Condition")

##### QC ####
human.t.2 <- qc(human.t)
human.u.2 <- qc(human.urt)
pig.t.2 <- qc(pig.t)
pig.u.2 <- qc(pig.urt)

VlnPlot(human.t.2, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.1) & theme(plot.title = element_text(size=10))

##### Step 3: Create Homo-gene only Human and Pig Matrices #####
human.t.homolog <- subset(human.t.2,
                          features = homolog_list)
human.urt.homolog <- subset(human.u.2,
                            features = homolog_list)
pig.t.homolog <- subset(pig.t.2,
                        features = homolog_list)
pig.urt.homolog <- subset(pig.u.2,
                          features = homolog_list)
### UMI MAPPING ###
# Treated #
# Find the unique and duplicated cells between the datasets
unique_human_cells <- setdiff(colnames(human.t.homolog@assays$RNA@counts), colnames(pig.t.homolog@assays$RNA@counts)) # 8 unique human cells
unique_porcine_cells <- setdiff(colnames(pig.t.homolog@assays$RNA@counts), colnames(human.t.homolog@assays$RNA@counts)) # 1767 unique porcine cells
duplicated_cells <- intersect(colnames(human.t.homolog@assays$RNA@counts), colnames(pig.t.homolog@assays$RNA@counts)) # 4224 duplicated cells

# Find the unique human cells
human.t.homolog.unique = human.t.homolog[, unique_human_cells]

# Merge the Seurat objects so duplicated cells are only included once (find the union)
merged.t.2 <- merge(human.t.homolog.unique, 
                  y = pig.t.homolog, 
                  project = "Merged")

# add metadata
merged.t.2@meta.data$origin <- NULL
merged.t.2@meta.data$percent_umi_human <- NULL
merged.t.2@meta.data$percent_umi_porcine <- NULL
merged.t.2@meta.data$origin[colnames(merged.t.2) %in% unique_human_cells] <- "human"
merged.t.2@meta.data$origin[colnames(merged.t.2) %in% unique_porcine_cells] <- "porcine"

# For unique human cells, set percent_umi_human to 100 and percent_umi_porcine to 0
merged.t.2@meta.data$percent_umi_human[colnames(merged.t.2) %in% unique_human_cells] <- 100
merged.t.2@meta.data$percent_umi_porcine[colnames(merged.t.2) %in% unique_human_cells] <- 0
merged.t.2@meta.data$percent_umi_porcine[colnames(merged.t.2) %in% unique_porcine_cells] <- 100
merged.t.2@meta.data$percent_umi_human[colnames(merged.t.2) %in% unique_porcine_cells] <- 0

# Calculate the total UMIs for each organism in each cell
human_umi_counts <- Matrix::colSums(human.t.homolog@assays$RNA@counts)
porcine_umi_counts <- Matrix::colSums(pig.t.homolog@assays$RNA@counts)

for (cell in duplicated_cells) {
  # fill in "origin" metadata with assigned organism
  if (human_umi_counts[cell] > porcine_umi_counts[cell]) {
    merged.t.2@meta.data$origin[colnames(merged.t.2@assays$RNA@counts) == cell] <- "human"
  } else {
    merged.t.2@meta.data$origin[colnames(merged.t.2@assays$RNA@counts) == cell] <- "porcine"
  }
    # fill in percent_umi_human with percentage of umi mapped to human
    total_umi <- human_umi_counts[cell] + porcine_umi_counts[cell]
    merged.t.2@meta.data$percent_umi_human[colnames(merged.t.2@assays$RNA@counts) == cell] <- round((human_umi_counts[cell] / total_umi) * 100, 2)
    # fill in percent_umi_human with percentage of umi mapped to porcine
    merged.t.2@meta.data$percent_umi_porcine[colnames(merged.t.2@assays$RNA@counts) == cell] <- round((porcine_umi_counts[cell] / total_umi) * 100, 2)
  }

table(merged.t.2@meta.data$origin)
# human porcine 
# 65    5934

# Look at UMI density
# convert the data to a data frame
df <- data.frame(
  percent_umi_human = merged.t.2@meta.data$percent_umi_human,
  percent_umi_porcine = merged.t.2@meta.data$percent_umi_porcine
)
# Add a new column 'color' to the dataframe based on the condition
df$color <- ifelse(df$percent_umi_porcine > 50, "red", "blue")

# Generate the plot
ggplot(df, aes(x = percent_umi_human, y = percent_umi_porcine, color = color)) +
  geom_point(alpha = 0.5) +
  scale_color_identity() +
  theme_minimal() +
  labs(
    x = "Percent UMI Human",
    y = "Percent UMI Porcine",
    title = "Treated Kidney 2 Percent UMI Human vs Percent UMI Porcine") + theme_bw()

# this means that 57 cells out of 1767 cells with porcine expression matrix needs to be substituted with human matrix

### create matrix with human and pig cells with assigned origin ###
human_cell_barcodes2 <- merged.t.2@assays$RNA@counts@Dimnames[[2]][merged.t.2@meta.data$origin == "human"]
human.cells2 <- human.t.homolog[, human_cell_barcodes2]
human.cells2@meta.data$origin <- rep("human", nrow(human.cells2@meta.data))

# Subset pig.t to exclude cells with matching barcodes
pig_cell_barcodes2 <- merged.t.2@assays$RNA@counts@Dimnames[[2]][merged.t.2@meta.data$origin == "porcine"]
pig.cells2 <- pig.t.homolog[, pig_cell_barcodes2]
pig.cells2@meta.data$origin <- rep("porcine", nrow(pig.cells2@meta.data))

treated2 <- merge(human.cells2, 
                 pig.cells2, 
                 project = "Merged")

# Untreated Room Temperature #
# Find the unique and duplicated cells between the datasets
unique_human_cells <- setdiff(colnames(human.urt.homolog@assays$RNA@counts), colnames(pig.urt.homolog@assays$RNA@counts)) # 0 cells
unique_porcine_cells <- setdiff(colnames(pig.urt.homolog@assays$RNA@counts), colnames(human.urt.homolog@assays$RNA@counts)) # 2624 cells
duplicated_cells <- intersect(colnames(human.urt.homolog@assays$RNA@counts), colnames(pig.urt.homolog@assays$RNA@counts)) # 5346 cells

# Find the unique human cells
human.urt.homolog.unique = human.urt.homolog[, unique_human_cells] # NONE

# Use pig.urt.homolog instead of merged object
merged.u.2 <- pig.urt.homolog
# add metadata
merged.u.2@meta.data$origin <- NULL
merged.u.2@meta.data$percent_umi_human <- NULL
merged.u.2@meta.data$percent_umi_porcine <- NULL
merged.u.2@meta.data$origin[colnames(merged.u.2) %in% unique_human_cells] <- "human"
merged.u.2@meta.data$origin[colnames(merged.u.2) %in% unique_porcine_cells] <- "porcine"

# For unique human cells, set percent_umi_human to 100 and percent_umi_porcine to 0
merged.u.2@meta.data$percent_umi_human[colnames(merged.u.2) %in% unique_human_cells] <- 100
merged.u.2@meta.data$percent_umi_porcine[colnames(merged.u.2) %in% unique_human_cells] <- 0
merged.u.2@meta.data$percent_umi_porcine[colnames(merged.u.2) %in% unique_porcine_cells] <- 100
merged.u.2@meta.data$percent_umi_human[colnames(merged.u.2) %in% unique_porcine_cells] <- 0

# Calculate the total UMIs for each organism in each cell
human_umi_counts <- Matrix::colSums(human.urt.homolog@assays$RNA@counts)
porcine_umi_counts <- Matrix::colSums(pig.urt.homolog@assays$RNA@counts)

for (cell in duplicated_cells) {
  # fill in "origin" metadata with assigned organism
  if (human_umi_counts[cell] > porcine_umi_counts[cell]) {
    merged.u.2@meta.data$origin[colnames(merged.u.2@assays$RNA@counts) == cell] <- "human"
  } else {
    merged.u.2@meta.data$origin[colnames(merged.u.2@assays$RNA@counts) == cell] <- "porcine"
  }
    # fill in percent_umi_human with percentage of umi mapped to human
    total_umi <- human_umi_counts[cell] + porcine_umi_counts[cell]
    merged.u.2@meta.data$percent_umi_human[colnames(merged.u.2@assays$RNA@counts) == cell] <- round((human_umi_counts[cell] / total_umi) * 100, 2)
    # fill in percent_umi_human with percentage of umi mapped to porcine
    merged.u.2@meta.data$percent_umi_porcine[colnames(merged.u.2@assays$RNA@counts) == cell] <- round((porcine_umi_counts[cell] / total_umi) * 100, 2)
  }

table(merged.u.2@meta.data$origin)
# porcine 
# 7970 

# Look at UMI density
# convert the data to a data frame
df <- data.frame(
  percent_umi_human = merged.u.2@meta.data$percent_umi_human,
  percent_umi_porcine = merged.u.2@meta.data$percent_umi_porcine
)
# Add a new column 'color' to the dataframe based on the condition
df$color <- ifelse(df$percent_umi_porcine > 50, "red", "blue")

# Generate the plot
ggplot(df, aes(x = percent_umi_human, y = percent_umi_porcine, color = color)) +
  geom_point(alpha = 0.5) +
  scale_color_identity() +
  theme_minimal() +
  labs(
    x = "Percent UMI Human",
    y = "Percent UMI Porcine",
    title = "Untreated Kidney 2 Percent UMI Human vs Percent UMI Porcine") + theme_bw() + xlim(0,100) + ylim(0,100)

untreated2 <- merged.u.2
#########################################################################
# final outputs:
# homolog-genes only:
# kidney 1: treated, untreated
# kidney 2: treated2, untreated2
# reference with full set of genes:
# kidney 1: human.t, human.u, pig.t, and pig.u
# kidney 2: human.t.2, human.u.2, pig.t.2, and pig.u.2
#########################################################################

#####################################
##### SCTransform + Integration #####
#####################################
library(glmGamPoi)
library(sctransform)
# perform normalization and dimensionality reduction
# normalize and run dimensionality reduction on control dataset
treated <- SCTransform(treated, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

untreated <- SCTransform(untreated, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

treated2 <- SCTransform(treated2, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

untreated2 <- SCTransform(untreated2, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

# Perform integration using pearson residuals
condition.list <- list(treated, untreated, treated2, untreated2)
features <- SelectIntegrationFeatures(object.list = condition.list, nfeatures = 3000)
condition.list <- PrepSCTIntegration(object.list = condition.list, anchor.features = features)

# integrate the two datasets
anchors <- FindIntegrationAnchors(object.list = condition.list, 
                                  normalization.method = "SCT",
                                  anchor.features = features)
kidneys <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Perform an integrated analysis
kidneys <- RunPCA(kidneys, verbose = FALSE)
kidneys <- RunUMAP(kidneys, reduction = "pca", dims = 1:30, verbose = FALSE)
kidneys <- FindNeighbors(kidneys, reduction = "pca", dims = 1:30)
kidneys <- FindClusters(kidneys, resolution = 0.3)

#### umap
DimPlot(kidneys, reduction = "umap")
# view by origin
DimPlot(kidneys, reduction = "umap", label = T, split.by = "origin")

#### cell type annotation ####
# Figure 1 D
DotPlot(kidneys, 
          features = c(
                       "SLC34A1", # proximal tubules cells,
                       "SLC6A18",
                       "SLC47A2",
                       
                       "STMN1", # proximal tubule VIM+
                       "PCLAF",
                       "HMGB2",
                       
                       "CALD1", # proximal tubule VIM+
                        "VIM",
                       "ANXA2",
                       
                       "SLC12A1", # thick ascending limb
                       "CLDN19",
                       #"UMOD",
                       "KRT7",
                       
                       "CALB1", # distal tubular cells
                       "AVPR2",
                       "AQP2",
                       
                       "SLC26A7", # IC type A
                       "SLC4A1",
                       "ATP6V0A4",
                       "HMX2", # IC type B
                       "LUM",
                       "SLC26A4",
                       "ATP6V1G3", # intercalated cells
                       "EMCN", # endo
                       "PECAM1",
                       "CDH5",
                       "HLA-DRA",
                       "LYZ",# Macrophage
                       "C1QA"
                       )) + theme(axis.text.x = element_text(angle = 90))
