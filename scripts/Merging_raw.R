## Data Import & Organization
library(here)             
library(data.table)
library(arrow)

## Data Manipulation
library(dplyr)            
library(tidyverse)        
## Spatial & Single-Cell Processing
library(Seurat) #latest version - v5          
# library(SeuratExtend)     #could'nt locate package
library(spacexr)          
library(glmGamPoi)        

## Integration & Batch Correction
library(harmony)          

## Visualization
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


#######SECTION 2: LOAD IN RDS OBJECTS, PERFORM QUALITY CONTROL, NORMALIZE DATA, DIMENTIONALITY REDUCTION AND CLUSTERING###########
# ------------------------------------------------------------------
# Important Note on Object Naming
# ------------------------------------------------------------------
# When loading in your saved .rds files, avoid using the same
# name as the FOV identifiers (the names you assigned to each field of view in the
# Xenium Analyzer) or orig.ident values from Xenium. Using distinct names helps
# prevent confusion and accidental overwriting.

# ------------------------------------------------------------------
# Example: Setting Working Directory and Reading an RDS File
# ------------------------------------------------------------------
# Set the working directory to the folder containing your saved RDS objects
setwd("D:/work/Xenium/Input_RDS") 


#############SECTION: MERGING DATA##################################################

# ------------------------------------------------------------------
# Multi-Sample Merge Example (Using raw RDS, unprocessed data)
# ------------------------------------------------------------------
#Read in original rds files (do not merge manipulated samples from section 2)

MAID <-readRDS("M_70_Seurat_obj.rds")
MAID <- UpdateSeuratObject(MAID)
TBI <-readRDS("M_38_Seurat_obj.rds")
TBI <- UpdateSeuratObject(TBI)
SEPSIS <-readRDS("M_56_Seurat_obj.rds")
SEPSIS <- UpdateSeuratObject(SEPSIS)
HIBI <-readRDS("M_44_Seurat_obj.rds")
HIBI <- UpdateSeuratObject(HIBI)

#assign orig IDs (FOVs from Xenium run)
MAID$orig.ident <- "M_70"
TBI$orig.ident <- "M_38"
SEPSIS$orig.ident <- "M_56"
HIBI$orig.ident <- "M_44"

##Add metadata columns for your conditions (e.g. developmental stage, sex, disease state etc)
#example here is Sex, but you can add several columns as needed
# Initialize a new metadata column 'Sex' with NA values for all cells
#merged_obj$Sex <- NA

# Assign 'Female' to cells from sample F1 and 'Male' to cells from sample M1 and M2
# This works by selecting cells where orig.ident matches the sample ID
#merged_obj$Sex[merged_obj$orig.ident %in% c("M1", "M2")] <- "Male"
#merged_obj$Sex[merged_obj$orig.ident %in% c("F1"] <- "Female"

# Merge samples using add.cell.ids to prefix cell barcodes
merged_obj <- merge(
  x = MAID,
  y = list(TBI, SEPSIS, HIBI),
  add.cell.ids = c("MAID", "TBI", "SEPSIS", "HIBI"))  

# Verify orig.ident values in metadata (should reflect FOV names assigned in Xenium Analyzer)
# - unique(): shows which categories are present
unique(merged_obj$orig.ident)

#to check if cell counts per sample are ok
# - table(): gives counts of each category to ensure all groups are represented
table(merged_obj$orig.ident)

#to save merged obj as rds
saveRDS(merged_obj, "merged_raw_obj.rds")


#to recall obj later
merged_obj <- readRDS("merged_raw_obj.rds") #"merged_raw_obj.rds"

##Add metadata columns for your conditions (e.g. developmental stage, sex, disease state etc)
#example here is Sex, but you can add several columns as needed
# Initialize a new metadata column 'Sex' with NA values for all cells
#merged_obj$Sex <- NA

# Assign 'Female' to cells from sample F1 and 'Male' to cells from sample M1 and M2
# This works by selecting cells where orig.ident matches the sample ID
#merged_obj$Sex[merged_obj$orig.ident %in% c("M1", "M2")] <- "Male"

# Simple check of 'Sex' metadata column values
# - unique(): shows which categories are present (expect "Male" and "Female")
# - table(): gives counts of each category to ensure both groups are represented
#unique(merged_obj$Sex)
# Count each category
#table(merged_obj$Sex)


#my metadata (for the first run) is condition type (TBI, MAID, Sepsis, HIBI), age at death
merged_obj$Condition <- NA
merged_obj$DeathAge <- NA

# Assign 'condition or age' to cells from respective cases 
# This works by selecting cells where orig.ident matches the sample IDs
#orig IDs are changed to FOV, so use that names
merged_obj$DeathAge[merged_obj$orig.ident %in% ("M_70")] <- "37"
merged_obj$DeathAge[merged_obj$orig.ident %in% ("M_38")] <- "58"
merged_obj$DeathAge[merged_obj$orig.ident %in% ("M_56")] <- "68"
merged_obj$DeathAge[merged_obj$orig.ident %in% ("M_44")] <- "74"

merged_obj$Condition[merged_obj$orig.ident %in% ("M_70")] <- "MAID"
merged_obj$Condition[merged_obj$orig.ident %in% ("M_38")] <- "TBI"
merged_obj$Condition[merged_obj$orig.ident %in% ("M_56")] <- "Sepsis"
merged_obj$Condition[merged_obj$orig.ident %in% ("M_44")] <- "HIBI"

# Simple check of 'Condition' metadata column values
# - unique(): shows which categories are present 
# - table(): gives counts of each category to ensure both groups are represented
unique(merged_obj$Condition)
# Count each category
table(merged_obj$Condition)

unique(merged_obj$DeathAge)
# Count each category
table(merged_obj$DeathAge)

#to save merged obj that contains metadata as rds
saveRDS(merged_obj, "merged_raw_metadata_obj.rds")

#Congrats you have now created your merged object with all your samples!
#And you know the workflow - follow steps 1-6 in Section 2 (QC, PCA, dimensionality reduction, clustering, visualization)

# ------------------------------------------------------------------
#QC, SCTransform, PCA, UMAP, Clustering - same as Single-Sample Analysis Workflow but done on the merged file
# ------------------------------------------------------------------
#

#to recall merged file 
merged_obj <- readRDS("merged_raw_metadata_obj.rds")

# 1. Quality Control 
# Visualize QC metrics
# - nFeature_Xenium: number of unique transcript features detected per cell (i.e., how many genes/transcripts)
# - nCount_Xenium: total transcript count per cell (sum of all detected transcript counts)

#quality control - feature means total genes picked which can be less than 
#the total panel genes 266 predesigned (human v1) + 100 custom (Kraus) 
# as some genes might not be expressed in our sample. 
#Count refer to number of transcripts for all the genes.

VlnPlot(merged_obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size =0)

# Filter cells based on QC thresholds (only if subsetting is needed e.g. for a region or transcripts of a specific size)
merged_obj <- subset(merged_obj, subset = nFeature_Xenium > 5 & nCount_Xenium > 0) 

#make Vln plot again
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
# Run the automated test
#dim_reso_test(
#  dataset = merged_obj,
#  dims_list = list(1:20, 1:30, 1:40, 1:50),  # test PC ranges
#  resolutions = c(0.1, 0.4, 0.8),    # test resolutions
#  # "Z:/Wellington Lab/Mehwish/Xenium_human_run1/B_HumanXenium_Run1_OG_rds"
#  output_dir = "D:/work/Xenium/Output", # where to save PDFs
#  red = "pca",             # dimensionality reduction to use
#  group_vars = c("orig.ident", "Condition")  # metadata for coloring
#)

#Choose best PC/resolution combo from PDFs and then rerun on main object
merged_obj <- FindNeighbors(merged_obj, dims = 1:30, reduction = "pca")
merged_obj <- FindClusters(merged_obj, resolution = 0.3)
merged_obj <- RunUMAP(merged_obj, dims = 1:30, reduction = "pca")

"GFAP" %in% rownames(merged_obj)
FeaturePlot(merged_obj, features = "SLC17A6")

# Option 2: PCA Heatmap (can also create these heatmaps after selecting the PCA above)
# - Use DimHeatmap() to visualize the top features (genes) driving each PC.
# - Adjust 'dims' to specify which PCs to examine; 'cells' controls how many cells to display per PC.
# - 'balanced = TRUE' scales positive and negative loadings equally.
DimHeatmap(merged_obj, dims = 1:30, cells = 500, balanced = TRUE)



# 4. Non-linear Embedding: UMAP
# UMAP parameters:
# - dims: PCs to include (e.g., 1:20)
merged_obj <- RunUMAP(merged_obj, dims = 1:30)

# 5. Clustering
# Identify clusters
# - resolution: higher => more clusters
merged_obj <- FindNeighbors(merged_obj, dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 0.3)

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
  fov = "M_44",
  cols = merged_obj@misc$cluster_colors,
  size = 0.5,
  border.size = NA,
  axes = TRUE,
  dark.background = TRUE
) + DarkTheme()



