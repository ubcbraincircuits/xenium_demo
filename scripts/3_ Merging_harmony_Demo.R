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


#Start by using your original merged object from section 3 
#ensure the merged obj did not go through any manipulation
merged_obj <- readRDS("merged_raw_obj.rds") #"merged_raw_obj.rds"


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
  group.by = "Condition"
) + DarkTheme() + coord_fixed() + ggtitle("merged_obj")

# Plot spatial FOV view
ImageDimPlot(
  merged_obj,
  fov = "X3_fov",
  cols = merged_obj@misc$cluster_colors,
  size = 0.5,
  border.size = NA,
  axes = TRUE,
  dark.background = TRUE
) + DarkTheme()

#When subseting integrated data, make sure switch assay to 'Xenium' and 
#follow block of code above: run SCT, PCA, integratelayers, dimensionality reduction, clustering etc. 