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
#Read in original rds files (do not merge manipulated samples from section 2)
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


# Merge samples using add.cell.ids to prefix cell barcodes
merged_obj <- merge(
  x = Mouse_2L_obj,
  y = list(Mouse_3L_obj, Mouse_4L_obj, Mouse_5L_obj),
  add.cell.ids = c("Region_2_left", "Region_3_left", "Region_4_left", "Region_5_left"))  

#to save merged obj as rds
saveRDS(merged_obj, "merged_raw_obj.rds")


#to recall obj later
merged_obj <- readRDS("merged_raw_obj.rds") #"merged_raw_obj.rds"


##Add metadata columns for your conditions (e.g. developmental stage, sex, disease state etc)
#example here is Orientation, but you can add several columns as needed
# Initialize a new metadata column 'Orientation' with NA values for all cells
merged_obj$Orientation <- NA

# Assign 'Left' to cells from sample X2_fov, X3_fov, X4_fov, X5_fov and 'Right' to cells from sample ...
# Note samples from Right region are not included here with regard to resources
# This works by selecting cells where orig.ident matches the sample ID
merged_obj$Orientation[merged_obj$orig.ident %in% c("X2_fov", "X3_fov")] <- "Left"
#merged_obj$Orientation[merged_obj$orig.ident %in% c("X7_fov"] <- "Right"

# Verify orig.ident values in metadata (should reflect FOV names assigned in Xenium Analyzer)
# - unique(): shows which categories are present
unique(merged_obj$orig.ident)
unique(merged_obj$Orientation)

#to check if cell counts per sample are ok
# - table(): gives counts of each category to ensure all groups are represented
table(merged_obj$orig.ident)
table(merged_obj$Orientation)


#to save merged obj that contains metadata as rds
saveRDS(merged_obj, "merged_raw_metadata_obj.rds")

#Congrats you have now created your merged object with all your samples!
#And you know the workflow - follow steps 1-6 in Section 2 (QC, PCA, dimensionality reduction, clustering, visualization)

# ------------------------------------------------------------------
#QC, SCTransform, PCA, UMAP, Clustering - same as Single-Sample Analysis Workflow but done on the merged file
# ------------------------------------------------------------------
#
# Clear memory from previous objects
rm(list = ls())
gc()

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

# Option 2: PCA Heatmap (can also create these heatmaps after selecting the PCA above)
# - Use DimHeatmap() to visualize the top features (genes) driving each PC.
# - Adjust 'dims' to specify which PCs to examine; 'cells' controls how many cells to display per PC.
# - 'balanced = TRUE' scales positive and negative loadings equally.
DimHeatmap(merged_obj, dims = 1:12, cells = 50, balanced = TRUE)


#Larissa's loop for testing resolution and dims for best selection
dim_reso_test <- function(dataset, dims_list, resolutions, output_dir, red, group_vars) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)  
  }
  for (dims in dims_list) {
    max_dim <- max(dims)
    for (res in resolutions) {
      dat_mod <- dataset 
      dat_mod <- FindNeighbors(dat_mod, dims = dims, reduction = red)
      dat_mod <- FindClusters(dat_mod, resolution = res)
      dat_mod <- RunUMAP(dat_mod, dims = dims, reduction = red)
      #assign colours
      current_clusters <- levels(Idents(dat_mod))
      # Generate distinct color palette (Polychrome handles up to 40–50 unique colors)
      if (length(current_clusters) <= 36) {
        palette <- Polychrome::palette36.colors(length(current_clusters))
      } else {
        # For >36 clusters, extend smoothly with hue palette
        palette <- scales::hue_pal()(length(current_clusters))
      }
      plot_pal <- palette[1:length(current_clusters)]
      names(plot_pal) <- current_clusters
      
      plots <- list(DimPlot(dat_mod, shuffle = TRUE, label = TRUE, cols = plot_pal) + coord_fixed())
      for (group_var in group_vars) {
        plots <- append(plots, list(DimPlot(dat_mod, group.by = group_var, shuffle = TRUE) + coord_fixed()))
      }
      # Create the vertical stack
      combined_plot <- patchwork::wrap_plots(plots, ncol = 1)
      
      # Add the metadata title
      plot_title <- paste("Dimensions:", max_dim, "| Resolution:", res)
      combined_plot <- combined_plot + patchwork::plot_annotation(title = plot_title)
      
      # Define PNG file path
      file_name <- paste0("Dim_", max_dim, "_Res_", res, ".png")
      png_path <- file.path(output_dir, file_name)
      
      # Save the plot
      # Height is dynamic based on number of plots; 300 DPI ensures crisp labels
      ggplot2::ggsave(
        filename = png_path,
        plot = combined_plot,
        width = 8,
        height = 5 * length(plots),
        dpi = 300,
        limitsize = FALSE
      )
    }
  }
}


# Run the automated test at desired dimension and resolution settings and examine results
dim_reso_test(
  dataset = merged_obj,
  dims_list = list(1:20, 1:30, 1:40),  # test PC ranges
  resolutions = c(0.1, 0.4),    # test resolutions
  # "Z:/Wellington Lab/Mehwish/Xenium_human_run1/B_HumanXenium_Run1_OG_rds"
  output_dir = "D:/work/Xenium/Output", # where to save PDFs
  red = "pca",             # dimensionality reduction to use
  group_vars = c("orig.ident", "Orientation")  # metadata for coloring
)


#Choose best PC/resolution combo from PDFs and then rerun on main object


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

# --- Plot UMAPs with locked colors ---
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

# --- Plot spatial FOV view with locked colors ---
ImageDimPlot(
  merged_obj,
  fov = "X3_fov",
  cols = merged_obj@misc$cluster_colors,
  size = 0.5,
  border.size = NA,
  axes = TRUE,
  dark.background = TRUE
) + DarkTheme()



