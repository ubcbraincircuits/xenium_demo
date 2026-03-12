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


####### LOAD IN RDS OBJECTS, PERFORM QUALITY CONTROL, NORMALIZE DATA, DIMENTIONALITY REDUCTION AND CLUSTERING ###########
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
setwd("D:/work/Github_demo/xenium_demo/Data")

# Load in the RDS object
Mouse_obj<-readRDS("Region_2_right-obj.rds")
# Update Seurat Objects to new structure for storing data/calculations. 
# For Seurat v3 objects, will validate object structure ensuring all keys 
# and feature names are formed properly.
Mouse_obj <- UpdateSeuratObject(Mouse_obj)


# Start by analyzing one sample at a time to understand the workflow and your data

# ------------------------------------------------------------------
# Single-Sample Analysis Workflow (QC, SCTransform, PCA, UMAP, Clustering)
# ------------------------------------------------------------------




# check ident and fov
table(Mouse_obj$orig.ident)
unique(Mouse_obj$orig.ident)

#Declare information like condition/Age or other categories before assigning them
Mouse_obj$Condition <- "Unknown"
#Mouse_obj$Condition[Mouse_obj$orig.ident %in% ("M_70")] <- "MAID"



# Metadata information is also included under:
# Mouse_obj@misc#run_metadata
# Run the above to see, or run View(Mouse_obj) to navigate the object.


## Preliminary analysis

# 1. Quality Control
# Visualize QC metrics

# Histogram 
# - nFeature_Xenium: number of unique transcript features detected per cell (i.e., how many genes/transcripts)
ggplot(Mouse_obj@meta.data, aes(x = nFeature_Xenium)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_classic() +
  #scale_y_log10() +        # for log scaling
  labs( #labels*
    title = "nFeature_Xenium",
    x = "Number of detected features",
    y = "Cell count"
  )

# - nCount_Xenium: total transcript count per cell (sum of all detected transcript counts)
# nCount_Xenium histogram
ggplot(Mouse_obj@meta.data, aes(x = nCount_Xenium)) +
  geom_histogram(bins = 25, fill = "firebrick", color = "black") +
  theme_classic() +
  #scale_y_log10() +        # for log scaling
  labs(
    title = "nCount_Xenium",
    x = "Total transcript counts",
    y = "Cell count"
  )


VlnPlot(Mouse_obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size =0)



# Filter cells based on QC thresholds
#Mouse_obj <- subset(Mouse_obj, subset = nFeature_Xenium > 2 & nCount_Xenium > 0) 
# Note: to compare multiple QC filters for downstream analysis
# assign filtered results to new object names 
# (e.g., obj1filtered) instead of overwriting obj1, obj2 etc.
#Mouse_obj_filtered <- subset(Mouse_obj, subset = nFeature_Xenium > 0 & nCount_Xenium > 0)


# View statistics of the object 
#n_genes <- dim(Mouse_obj)[1]
#n_cells <- dim(Mouse_obj)[2]


# 2. Normalization and Feature Selection with SCTransform
# SCTransform performs normalization, variance stabilization, and identifies variable features
Mouse_obj <- SCTransform(Mouse_obj, assay = "Xenium")

#You can choose the normalization method and other details this way:
# 02_run_seurat.R lines 11-109
# Link to documentation on SCTransform: https://satijalab.org/seurat/articles/sctransform_vignette.html
# Link to paper on SCTransform: https://link.springer.com/article/10.1186/s13059-021-02584-9

# function to plot residual variance vs gene expression 
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
    xlab("Gene mean") +
    ylab("Residual variance")
  if (annotate) {
    p <- p + geom_text_repel(
      data = subset(gene_var, gene %in% topn), aes(label = gene), color = "gray25",
      size = 1.8,
      nudge_y = 230 - subset(gene_var, gene %in% topn)[, col],
      direction = "x",
      angle = 90,
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

# test normalization by plotting "residual variance" vs "gene expression" 
# SCTResults is a way to pull data from SCTAssay object:
# documentation: https://satijalab.org/seurat/reference/sctresults#:~:text=Arguments,just%20returns%20the%20slot%20directly).

gene_attr_maid <- SCTResults(Mouse_obj, slot = "feature.attributes", assay = "SCT")

residualVarPlot(gene_attr_maid, max_resvar = 10, pt_size = 4)



# 3. Dimensionality Reduction: PCA
# PCA reduces high-dimensional data to principal components
# - features: use VariableFeatures from SCTransform
# - npcs: number of PCs to compute (e.g., 50)
# - Balance is Key: 
# - Including more PCs captures more variation in the data, but may introduce noise.
# - Including fewer PCs reduces noise but may risk missing important biological signals.

Mouse_obj <- RunPCA(Mouse_obj, features = VariableFeatures(Mouse_obj), npcs = 50, verbose = FALSE)


# Option 1: Examine ElbowPlot to choose how many PCs to use downstream

# In this dataset, the elbow appears around PCs 12-15. 
# This suggests that most of the biologically relevant variation is captured within the first set of components.
# PC's after this point may be mostly noise. Lower PC dimensions may decrease this noise but could also miss
# biologically relevant information.
# PCA is about maximizing variance in the first components. However, 
# variance is a quantity that makes most sense for continuous and non-sparse data.
# if you dont see a clear elbow, you data might be unsuitable for an elbow plot examination
ElbowPlot(Mouse_obj, ndims = 50)

# Option 2: PCA Heatmap
# - Use DimHeatmap() to visualize the top features (genes) driving each PC.
# - Adjust 'dims' to specify which PCs to examine; 'cells' controls how many cells to display per PC.
# - 'balanced = TRUE' scales positive and negative loadings equally.
DimHeatmap(Mouse_obj, dims = 1:12, cells = 500, balanced = TRUE)



# 4. Non-linear Embedding: UMAP
# UMAP parameters:
# - dims: PCs to include (e.g., 1:20)

Mouse_obj <- FindNeighbors(Mouse_obj, dims = 1:24)
Mouse_obj <- FindClusters(Mouse_obj, resolution = 0.1)
Mouse_obj <- RunUMAP(Mouse_obj, dims = 1:24)


## Set single colour pallete
# Get all cluster IDs
global_clusters <- levels(Idents(Mouse_obj))
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
Mouse_obj@misc$cluster_colors <- cluster_colors

# Ensure Idents is a factor with consistent levels
Idents(Mouse_obj) <- factor(Idents(Mouse_obj), levels = global_clusters)

#Plot UMAP (recommend adding PC and resol into the title to keep track of diff umaps)
DimPlot(Mouse_obj, group.by = "seurat_clusters", 
        reduction = "umap",
        cols = Mouse_obj@misc$cluster_colors) + DarkTheme() + coord_fixed() + ggtitle("UMAP REDUCTION")


#Plot your UMAP by metadata (ex. orig.ident, sex, stage, disease condition etc.)
#this is only after you have merged and/or integrated your data and added metadata columns (Section 3 and 4)
#DimPlot(Mouse_obj, reduction = "umap", cols = Mouse_obj@misc$cluster_colors) + DarkTheme() + coord_fixed() + ggtitle("UMAP REDUCTION grouped"



# 5. Clustering
# Identify clusters

# - resolution: higher => more clusters
Mouse_Neighbors <- FindNeighbors(Mouse_obj, dims = 1:20)
Mouse_Clustered <- FindClusters(Mouse_Neighbors, resolution = 0.1)


#ImageDimPlot(Mouse_Clustered, fov = "M1", cols = Mouse_obj@misc$cluster_colors, size = 0.5, border.size = NA,
#axes = TRUE, dark.background = TRUE) + DarkTheme()
ImageDimPlot(Mouse_Clustered, molecules = "Slc17a7", nmols = 10000, alpha = 0.3, mols.cols = "red")

# Note: After running SCTransform, it sets the active assay to "SCT" (normalized data).
# For FeaturePlots() and VlnPlots() plot raw Xenium counts, temporarily switch the default assay.
DefaultAssay(Mouse_Clustered) <- "Xenium"
#FeaturePlot: gene expression in UMAP space 
FeaturePlot(Mouse_obj, features = c("Aqp4"), label = TRUE)

#Violin Plot of raw counts (optional: set log = TRUE for log scale)
VlnPlot(Mouse_obj, features = c("Slc17a7", "Gfap", "Sla"), pt.size = 0, log = TRUE)
#uses the polychrome colours
VlnPlot_scCustom(Mouse_obj, features = c("Slc17a7", "Gfap", "Sla"), pt.size = 0, log = TRUE)

DotPlot(object = Mouse_obj, features = c("Aqp4", "Paqr5", "Trem2"), dot.min  = 0.1,
        dot.scale= 6, group.by = "nFeature_Xenium") +
  scale_color_gradientn(colors = c("steelblue", "white", "firebrick")) +
  theme_classic() + ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))
ImageFeaturePlot(Mouse_obj, fov = "X8fov", features = c("Aqp4", "Paqr5", "Trem2"))
#for single gene
ImageFeaturePlot(Mouse_obj, fov = "X8fov", features = ("Slc17a7"))





# ------------------------------------------------------------------
# Subsetting Data Based on Clusters of Interest
# ------------------------------------------------------------------
# If you are aiming to analyze a specific cell type or region you most likely need to subset your data
# Use a combination of marker gene expression and spatial location
# to identify the clusters you want to analyze further.
# you may need to perform several subsets to get to your cell type or region of interest based on how large your dataset is

###Several ways to subset (dependent on dataset)
##Option 1: Subset the object by clusters of interest (most common):
obj_sub <- subset(Mouse_obj, idents = c(0, 2, 4, 6))

##Option 2: subset based on marker expression thresholds (e.g. cells that express a marker gene or combination of them)
# Example: keep cells where MarkerGene1 expression > 1 in the subset
#obj_sub <- subset(Mouse_obj, subset = SLC17A7 > 1)

##Option 3: Crop your image down to a region of interest 
#Copy your original object
#obj1_crop <- Mouse_obj

#Crop the FOV (“M1”) by specifying x/y ranges
#(replace the numbers with your desired coordinates)
#every sample will have a diff set of coordinates so you use ImageDimPlot() for x and y coordinates
#in ImageDimPlot() what is labeled as y on the axis is actually x coordinates in the code,
#and what is labeled as x on the axis is actually y coordinates
#obj1_crop[["M_70"]] <- Crop(
#  obj1_crop[["M_70"]], x = c(100, 2000), y = c(500, 3000))

# Pull out all cell barcodes within that cropped region
#keep_cells <- Cells(obj1_crop[["M_70"]])

#Subset so that only those cropped cells remain
#obj_sub <- subset(obj1_crop, cells = keep_cells)

#Quick spatial check to confirm cropping
ImageDimPlot(obj_sub, fov = "X8fov", cols = Mouse_obj@misc$cluster_colors, size = 0.5, border.size = NA,
             axes = TRUE, dark.background = TRUE) + DarkTheme() 
ImageDimPlot(Mouse_obj, fov = "X8fov", cols = Mouse_obj@misc$cluster_colors, size = 0.5, border.size = NA,
             axes = TRUE, dark.background = TRUE) + DarkTheme() 

##Once you have your subset using  any one of the options above, you can perform steps 1-7 as above
# Before re-running SCTransform on the refined subset, switch default assay
DefaultAssay(obj_sub) <- "Xenium"
