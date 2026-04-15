## Data Import & Organization
library(here)             
library(data.table)
library(arrow)

## Data Manipulation
library(dplyr)            
library(tidyverse)        

## Spatial & Single-Cell Processing
library(Seurat) #latest version - v5          
library(spacexr)          
library(glmGamPoi)        

## Integration & Batch Correction
library(harmony)          

## Visualization
library(ggrepel)
library(ggplot2)          
library(patchwork)        
library(RColorBrewer)     
library(scCustomize)   
library(scattermore)
library(paletteer)    
library(Polychrome)

## Parallelization & Performance
library(future)  

#SingleR Analysis
library(SingleR)
library(rlang)
library(geometry)
library(purrr)

# Increase maximum future global size for handling large objects
options(future.globals.maxSize = 250000 * 1024^2)

####### LOAD IN RDS OBJECTS

# Set the working directory to the folder containing your saved RDS objects
setwd("D:/work/Github_demo/xenium_demo/Data")

############# SECTION: MERGING DATA##################################################

# ------------------------------------------------------------------
# Multi-Sample Merge Example (Using raw RDS, unprocessed data)
# ------------------------------------------------------------------
#Read in the rds files 
# Note: memory usage scales with the number of objects loaded. If you are running
#       less than 16GB of system memory, we suggest commenting out object 4 and 5.
Mouse_2L_obj <-readRDS("Region_2_left-obj.rds")
Mouse_2L_obj <- UpdateSeuratObject(Mouse_2L_obj)
Mouse_3L_obj <-readRDS("Region_3_left-obj.rds")
Mouse_3L_obj <- UpdateSeuratObject(Mouse_3L_obj)
Mouse_4L_obj <-readRDS("Region_4_left-obj.rds")
Mouse_4L_obj <- UpdateSeuratObject(Mouse_4L_obj)
Mouse_5L_obj <-readRDS("Region_5_left-obj.rds")
Mouse_5L_obj <- UpdateSeuratObject(Mouse_5L_obj)

# set fov's
# By manually creating an orig.ident column in the metadata and filling it with the FOV name you create
# a text label that links every individual cell's data profile to the physical tissue slice it came from
# which will help us differentiate them in downstream analysis
Mouse_2L_obj$orig.ident <- "X2_fov"
Mouse_3L_obj$orig.ident <- "X3_fov"
Mouse_4L_obj$orig.ident <- "X4_fov"
Mouse_5L_obj$orig.ident <- "X5_fov"


# We recommend running the following command, to familiarize yourself with a 
# single Seurat object
View(Mouse_2L_obj)

# Merge samples using add.cell.ids to prefix cell barcodes
merged_obj <- merge(
  x = Mouse_2L_obj,
  y = list(Mouse_3L_obj, Mouse_4L_obj, Mouse_5L_obj),
  add.cell.ids = c("Region_2_left", "Region_3_left", "Region_4_left", "Region_5_left"))  

# check fov names:
Images(merged_obj) # or use names(Mouse_obj@images)

# Contrast it to the merged object:
View(merged_obj)

# Save merged raw obj as rds to re-use later down the line
# Whenever we want to refresh or do new analysis.

# It is saved in the working directory
saveRDS(merged_obj, "merged_raw_obj.rds")

# We can also clear old objects
# Remove the objects from the environment
rm(Mouse_2L_obj, Mouse_3L_obj, Mouse_4L_obj, Mouse_5L_obj)
gc()

#to recall obj later
merged_obj <- readRDS("merged_raw_obj.rds") #"merged_raw_obj.rds"



##Add metadata columns for your conditions (e.g. developmental stage, sex, disease state etc)
# example here is Group, but it is only a place-holder for the purposes of 
# demonstration. It is not necessarily descriptive of the data. You can change
# add other fields, but atleast one field will be required
# for some visualizations downstream

# Initialize a new metadata column 'Group' with NA values for all cells
merged_obj$Group <- NA

# Assign 'A' to cells from sample X2_fov, X3_fov, and 'B' to cells from sample X4_fov, X5_fov 
# This works by selecting cells where orig.ident matches the sample ID
merged_obj$Group[merged_obj$orig.ident %in% c("X2_fov", "X3_fov")] <- "A"
merged_obj$Group[merged_obj$orig.ident %in% c("X4_fov", "X5_fov")] <- "B"

# Verify orig.ident values in metadata (should reflect FOV names assigned in Xenium Analyzer)
# - unique(): shows which categories are present
unique(merged_obj$orig.ident)
unique(merged_obj$Group)

#to check if cell counts per sample are ok
# - table(): gives counts of each category to ensure all groups are represented
table(merged_obj$orig.ident)
table(merged_obj$Group)

#to save merged obj that contains metadata as rds
saveRDS(merged_obj, "merged_raw_metadata_obj.rds")

#Congrats you have now created your merged object with all your samples!
#And you know the workflow - follow steps 1-6 in Section 2 (QC, PCA, dimensionality reduction, clustering, visualization)




# ------------------------------------------------------------------
#QC, SCTransform, PCA, UMAP, Clustering - same as Single-Sample Analysis Workflow but done on the merged file
# ------------------------------------------------------------------
#

# Recall merged file 
merged_obj <- readRDS("merged_raw_metadata_obj.rds")

# 1. Quality Control 
# Visualize QC metrics
VlnPlot(merged_obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size =0)

# Filter cells based on QC thresholds (only if subsetting is needed e.g. for a region or transcripts of a specific size)
merged_obj <- subset(merged_obj, subset = nFeature_Xenium > 5 & nCount_Xenium > 0) 

# make Vln plot again
VlnPlot(merged_obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size =0)



  # 2. Normalization and Feature Selection with SCTransform
# SCTransform performs normalization, variance stabilization, and identifies variable features
merged_obj <- SCTransform(merged_obj, assay = "Xenium")

# function to plot residual variance vs gene expression:  to visualize variance stabilization
# updated for merged obj
residualVarPlot <- function(gene_var, xaxis = "gmean", max_resvar = 100, ntop = 20, annotate = F, pt_size = 1.1) {
  gene_var$gene <- rownames(gene_var)
  topn <- subset(gene_var, rank(-gene_var[, "residual_variance"]) <= ntop)$gene
  gene_var[gene_var$residual_variance > max_resvar, "residual_variance"] <- max_resvar
  p <- ggplot(gene_var, aes_string(xaxis, "residual_variance")) +
    geom_scattermore(pointsize = pt_size, shape = 16, alpha = 0.5, color = "#43a2ca") +
    geom_scattermore(data = subset(gene_var, gene %in% topn), pointsize = pt_size, shape = 16, alpha = 1.0, color = "deeppink") +
    geom_hline(yintercept = 1, color = "#4daf4a", size = 0.9, linetype = "dashed") +
    geom_smooth(method = "loess", span = 0.1, size = 0.9, formula = "y ~ x", color = "#e41a1c") +
    scale_y_continuous(trans = "sqrt", breaks = c(0, 1, 10, 25, 50, 100, 150), limits = c(0, max_resvar + 1)) +
    scale_x_continuous(trans = "log10", breaks = c(0.001, 0.01, 0.1, 1, 10, 100), labels = MASS::rational) + #, 100
    # facet_wrap(~ model, ncol=3, scales = 'free_y') +
    xlab("Gene mean expression") +
    ylab("Residual variance")
  if (annotate) {
    p <- p + geom_text_repel(
      data = subset(gene_var, gene %in% topn), aes(label = gene), color = "gray25",
      size = 1.8,
      direction = "x",
      #angle = 90,
      vjust = 0.5,
      hjust = 0.5,
      segment.size = 0.2,
      segment.alpha = 0.2
    )
  }
  
  cat("Top", ntop, "genes by residual variance:\n")
  print(topn)
  
  return(p)
}

# Extracts and combines SCT attributes for merged objects
get_gene_attributes <- function(obj, assay = "SCT") {
  # Get the results (returns a list for merged objects)
  res <- SCTResults(obj, slot = "feature.attributes", assay = assay)
  
  # If it's already a single data frame, return it
  if (is.data.frame(res)) {
    res$gene <- rownames(res)
    res$model <- "Model_1"
    return(res)
  }
  
  # If it's a list, bind all data frames together
  combined_df <- do.call(rbind, lapply(names(res), function(model_name) {
    df <- res[[model_name]]
    df$gene <- rownames(df)
    df$model <- model_name # Keep track of which FOV/Sample the data comes from
    return(df)
  }))
  
  return(combined_df)
}

# test normalization by plotting "residual variance" vs "gene expression" 
# SCTResults is a way to pull data from SCTAssay object:
# documentation: https://satijalab.org/seurat/reference/sctresults#:~:text=Arguments,just%20returns%20the%20slot%20directly).
gene_attr_merged <- get_gene_attributes(merged_obj, assay = "SCT")

residualVarPlot(gene_attr_merged, max_resvar = 15, pt_size = 3, annotate = T)


# 3. Dimensionality Reduction: PCA
# PCA reduces high-dimensional data to principal components
# - features: use VariableFeatures from SCTransform
# - npcs: number of PCs to compute (e.g., 50)
# - Balance is Key: 
# - Including more PCs captures more variation in the data, but may introduce noise.
# - Including fewer PCs reduces noise but may risk missing important biological signals.

merged_obj <- RunPCA(merged_obj, features = VariableFeatures(merged_obj), npcs = 50, verbose = FALSE)

Reductions(merged_obj)

# Option 1: Examine ElbowPlot to choose how many PCs to use downstream
ElbowPlot(merged_obj, ndims = 50) 

# Alternatively, we can also Visualize this as Percentage Variance Explained, 
# instead of change in Standard deviation this way:

# Extract standard deviations and calculate variance (stdev^2)
pca_var <- merged_obj[["pca"]]@stdev^2

# Calculate total variance from the scaled data layer
total_var <- sum(rowVars(LayerData(merged_obj, assay = "SCT", layer = "scale.data")))

# Create a dataframe for the first 50 PCs
plot_data <- data.frame(PC = 1:50, Pct_Variance = (pca_var[1:50] / total_var) * 100)

# Generate the Elbow Plot for Variance Explained
ggplot(plot_data, aes(x = PC, y = Pct_Variance)) +
  geom_point(size = 2) +
  labs(
    title = "Elbow Plot: Variance Explained",
    x = "Principal Component", 
    y = "% Variance explained"
  ) +
  theme_classic()

# Option 2: PCA Heatmap (can also create these heatmaps after selecting the PCA above)
# - Use DimHeatmap() to visualize the top features (genes) driving each PC.
# - Adjust 'dims' to specify which PCs to examine; 'cells' controls how many cells to display per PC.
# - 'balanced = TRUE' scales positive and negative loadings equally.
DimHeatmap(merged_obj, dims = 1:12, cells = 50, balanced = TRUE)

###############
# Recommended follow Larissa's loop from previous script.
# We skip that in interest of time
###############



# 4. Clustering
# Identify clusters
# - resolution: higher => more clusters
merged_obj <- FindNeighbors(merged_obj, dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 0.3)



# 5. Non-linear Embedding: UMAP
# UMAP parameters:
# - dims: PCs to include (e.g., 1:20)
merged_obj <- RunUMAP(merged_obj, dims = 1:30)



# 6.Visualizations

# Standardising colors
# Get all cluster IDs
global_clusters <- levels(Idents(merged_obj))
n <- length(global_clusters)

# Generate distinct color palette (Polychrome handles up to 40–50 unique colors)
if (n <= 36) {
  palette <- Polychrome::palette36.colors(n)
} else {
  # For >36 clusters, extend smoothly with hue palette
  palette <- scales::hue_pal()(n)
}

# Create a named color vector for consistent mapping
cluster_colors <- setNames(palette, global_clusters)

# Store mapping for reproducibility
merged_obj@misc$cluster_colors <- cluster_colors

# Ensure Idents is a factor with consistent levels
Idents(merged_obj) <- factor(Idents(merged_obj), levels = global_clusters)

# Plot UMAPs
DimPlot(
  merged_obj,
  reduction = "umap",
  cols = merged_obj@misc$cluster_colors
) + DarkTheme() + coord_fixed() + ggtitle("merged_obj")

DimPlot(
  merged_obj,
  reduction = "umap",
  cols = merged_obj@misc$cluster_colors,
  group.by = "Group"
) + DarkTheme() + coord_fixed() + ggtitle("merged_obj")

# Plot spatial FOV view
ImageDimPlot(
  merged_obj,
  fov = "X3fov",
  cols = merged_obj@misc$cluster_colors,
  size = 0.5,
  border.size = NA,
  axes = TRUE,
  dark.background = TRUE
) + DarkTheme()


# You can also perform subseting as discussed in previous section!




############# SECTION: HARMONY INTEGRATION ####################################################

# --- Merge vs. Integration: When to Use Which ---
# Merge (Seurat::merge):
#   • Combines multiple Seurat objects into one.
#   • Ideal when samples share similar technical conditions and batch effects are minor.
# Integration (Harmony):
#   • Actively corrects technical or batch differences between datasets.
#   • Essential when combining data from different experimental dates, sample preps, timepoints.
#   • Aligns shared biological signals across batches, improving clustering and visualization.
#   • For some datasets, it can wash away TRUE biological differences (e.g. integration healthy vs disease conditions).
# Decide what is best for your dataset!!

# For more detailed coverage, paper on Harmony: https://www.nature.com/articles/s41592-019-0619-0

# Note that for this dataset, Harmony is not quite applicable for our data, as it has
# very similar samples, and low batch effects. The purpose of the following script
# is to showcase how Harmony integration would be performed

# First lets clear previously loaded items
rm(list = ls())
gc()


# And using original merged object from the previous section 
# ensure the merged obj did not go through any manipulation
merged_obj <- readRDS("merged_raw_metadata_obj.rds")

#Normalize using SCTransform 
int_obj<- SCTransform(merged_obj,assay = "Xenium", verbose = FALSE)

#Run PCA 
int_obj<- RunPCA(int_obj,features = rownames(merged_obj), npcs = 100,
                 verbose = FALSE)

#Select PCs via ElbowPlot (or other methods discussed above)
ElbowPlot(int_obj,ndims = 100) 

# Integrate samples with Harmony
#    - orig.reduction: PCA embeddings to correct
#    - new.reduction: stores the batch-corrected embeddings
#    - normalization.method: use SCT-normalized data
int_obj<- IntegrateLayers(
  object = int_obj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = TRUE)

#Run UMAP, find neighbours and clustering
int_obj<- RunUMAP(int_obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
int_obj<- FindNeighbors(int_obj, reduction = "harmony", dims = 1:20)
int_obj<- FindClusters(int_obj, resolution = 0.1)



#Perform visualization and subseting as discussed in previous sections!

# Standardising colors
# Get all cluster IDs
global_clusters <- levels(Idents(merged_obj))
n <- length(global_clusters)

# Generate distinct color palette (Polychrome handles up to 40–50 unique colors)
if (n <= 36) {
  palette <- Polychrome::palette36.colors(n)
} else {
  # For >36 clusters, extend smoothly with hue palette
  palette <- scales::hue_pal()(n)
}

# Create a named color vector for consistent mapping
cluster_colors <- setNames(palette, global_clusters)

# Store mapping for reproducibility
int_obj@misc$cluster_colors <- cluster_colors

# Ensure Idents is a factor with consistent levels
Idents(int_obj) <- factor(Idents(int_obj), levels = global_clusters)

# Plot UMAPs
DimPlot(
  int_obj,
  reduction = "umap",
  cols = int_obj@misc$cluster_colors
) + DarkTheme() + coord_fixed() + ggtitle("int_obj")

DimPlot(
  int_obj,
  reduction = "umap",
  cols = int_obj@misc$cluster_colors,
  group.by = "Group"
) + DarkTheme() + coord_fixed() + ggtitle("int_obj")

# Plot spatial FOV view
ImageDimPlot(
  int_obj,
  fov = "X3fov",
  cols = int_obj@misc$cluster_colors,
  size = 0.5,
  border.size = NA,
  axes = TRUE,
  dark.background = TRUE
) + DarkTheme()

#When subseting integrated data, make sure switch assay to 'Xenium' and 
#follow block of code above: run SCT, PCA, integratelayers, dimensionality reduction, clustering etc. 


