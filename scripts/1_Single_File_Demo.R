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


####### SECTION:  LOAD IN RDS OBJECTS, PERFORM QUALITY CONTROL, NORMALIZE DATA, DIMENTIONALITY REDUCTION AND CLUSTERING ###########
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
# should be set to your downloads folder, eg: 
# setwd("C:/User/User/Downloads/Data")
setwd("D:/work/Github_demo/xenium_demo/Data")

# Load in the RDS object
Mouse_obj<-readRDS("Region_2_right-obj.rds")
# Update Seurat Objects to new structure for storing data/calculations. 
# For Seurat v3 objects, will validate object structure ensuring all keys 
# and feature names are formed properly.
Mouse_obj <- UpdateSeuratObject(Mouse_obj)


# Start by analyzing one sample at a time to understand the workflow and your data

# We recommend running the following command, to familiarize yourself with a 
# Seurat object
View(Mouse_obj) 
# note: The View function is case sensitive 
# View() - (capital 'V') is the built-in RStudio function that opens an 
#           interactive spreadsheet-style data viewer in a new tab.
# view() - (lowercase 'v') is a function from the tidyverse/tibble package. It
#           shows a static table style output in a new tab.

# You can look under Assays. It shows the different assays in the object
# with the number of [columns x rows]


# ------------------------------------------------------------------
# Single-Sample Analysis Workflow (QC, SCTransform, PCA, UMAP, Clustering)
# ------------------------------------------------------------------


# In a standard single-cell RNA-seq Seurat object, you only have matrices of cells and genes.
# However, in a spatial Seurat object, you need a place to store physical coordinates.
# Seurat created the @images slot specifically to hold FOV objects.

# The fov object is a container for storing coordinates of spatially-resolved single cells.
# Capable of storing multiple cell segmentation boundary masks.

# Your Seurat object might have one FOV, or it might contain many of distinct FOVs if 
# you imaged multiple separate regions of interest (like different brain regions) on the same slide

# check fov names:
Images(Mouse_obj) # or use names(Mouse_obj@images)


# You can add additional fields to the seurat object including supplementary information such as metadata.

# To include this, you must first Declare field name like condition/Age or other categories
#Mouse_obj$Condition <- "Unknown"

# Then assign a value to that field
#Mouse_obj$Condition[Mouse_obj$orig.ident %in% ("x8fov")] <- "MAID"



# Metadata information for the given data is also already included under:
# Mouse_obj@misc#run_metadata
# Run the above to see, or run View(Mouse_obj) to navigate the object.


## Preliminary analysis

# 1. Quality Control
# Visualize QC metrics
# - nFeature_Xenium: number of unique transcript features detected per cell (i.e., how many genes/transcripts)
# - nCount_Xenium: total transcript count per cell (sum of all detected transcript counts)

# quality control - feature means total genes picked which can be less than 
# the total panel genes 266 predesigned (human v1) + 100 custom (Kraus) 
# as some genes might not be expressed in our sample. 
# Count refer to number of transcripts for all the genes.

## Note that the data available has already been filtered for nFeature_Xenium > 0 & nCount_Xenium > 0

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


# include comments on recycling to this part

# Filter cells based on QC thresholds
#Mouse_obj <- subset(Mouse_obj, subset = nFeature_Xenium > 2 & nCount_Xenium > 0) 
# Note: to compare multiple QC filters for downstream analysis
# assign filtered results to new object names 
# (e.g., obj1filtered) instead of overwriting obj1, obj2 etc.
#Mouse_obj_filtered <- subset(Mouse_obj, subset = nFeature_Xenium > 0 & nCount_Xenium > 0)


# View statistics of the object 
n_genes <- dim(Mouse_obj)[1]
n_cells <- dim(Mouse_obj)[2]


# 2. Normalization and Feature Selection with SCTransform
# SCTransform performs normalization, variance stabilization, and identifies variable features
# Variance stabilization aims to remove technical variation (such as sequencing depth differences) 
# while preserving biological variance. This makes the data more comparable across cells.
Mouse_obj <- SCTransform(Mouse_obj, assay = "Xenium")

#You can choose the normalization method and other details this way:
# 02_run_seurat.R lines 11-109
# Link to documentation on SCTransform: https://satijalab.org/seurat/articles/sctransform_vignette.html
# Link to paper on SCTransform: https://link.springer.com/article/10.1186/s13059-021-02584-9


# For the 2 below function definitions, It is recommended to click the small downwards
# arrow to the left on line (tentative) 190 and 224, to remove the clutter.
# you can expand the function definition if you want insight into how these work.

# function to plot residual variance vs gene expression:  to visualize variance stabilization
# - gene_var:   The data frame containing gene attributes. 
# - xaxis:      (Default: "gmean") The name of the column to plot on the x-axis.
#               (gmean: geometric mean of gene expression.) 
# - max_resvar: The maximum limit for the y-axis (residual variance). (Default: 100) 
# - ntop:       (Default: 20) The number of top most variable genes to highlight. 
# - annotate:   (Default: FALSE). If set to TRUE, 
#               the function will draw text labels (gene names) next to the 'ntop' highlighted genes.
# - pt_size:    (Default: 1.1) The size of the dots drawn on the scatter plot.
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
gene_attr_mouse <- get_gene_attributes(Mouse_obj, assay = "SCT")

residualVarPlot(gene_attr_mouse, max_resvar = 10, pt_size = 3, annotate = T)

# 3. Dimensionality Reduction: PCA
# PCA reduces high-dimensional data to principal components with highest explained variance 
# - features: use VariableFeatures from SCTransform
# - npcs: number of PCs to compute (e.g., 50)
# - Balance is Key: 
# - Including more PCs captures more variation in the data, but may introduce noise.
# - Including fewer PCs reduces noise but may risk missing important biological signals.

# For an intuitive explanation of PCA, see this StatQuest video: https://www.youtube.com/watch?v=FgakZw6K1QQ

Mouse_obj <- RunPCA(Mouse_obj, features = VariableFeatures(Mouse_obj), npcs = 50, verbose = FALSE)


# Option 1: Examine ElbowPlot to choose how many PCs to use downstream

# In this dataset, the elbow appears around PCs 15-25. 
# This suggests that most of the biologically relevant variation is captured within the first set of components.
# PC's after this point may be mostly noise. Lower PC dimensions may decrease this noise but could also miss
# biologically relevant information.
# PCA is about maximizing variance in the first components. However, 
# variance is a quantity that makes most sense for continuous and non-sparse data.
# if you dont see a clear elbow, you data might be unsuitable for an elbow plot examination
ElbowPlot(Mouse_obj, ndims = 50)

# Alternatively, we can also Visualize this as Percentage Variance Explained, 
# instead of change in Standard deviation this way:

# Extract standard deviations and calculate variance (stdev^2)
pca_var <- Mouse_obj[["pca"]]@stdev^2

# Calculate total variance from the scaled data layer
total_var <- sum(rowVars(LayerData(Mouse_obj, assay = "SCT", layer = "scale.data")))

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


# Option 2: PCA Heatmap
# - Use DimHeatmap() to visualize the top features (genes) driving each PC.
# - Adjust 'dims' to specify which PCs to examine; 'cells' controls how many cells to display per PC.
# - 'balanced = TRUE' scales positive and negative loadings equally.
DimHeatmap(Mouse_obj, dims = 1:12, cells = 500, balanced = TRUE)


# loop for testing resolution and dims for best selection (Credit: Larissa Kraus)
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
# This is optional, and meant to be a "set it and forget it" test to be exhaustive
# and run overnight/over a weekend to help you choose a good resolution and
# number of dimensions. You can skip this and see the results on the resources page
# on Github
dim_reso_test(
  dataset = Mouse_obj,
  dims_list = list(1:20, 1:30, 1:40, 1:50),  # test PC ranges
  resolutions = c(0.1, 0.3, 0.5, 0.7),    # test resolutions
  output_dir = "D:/work/Github_demo/xenium_demo/Data/Output", # where to save PDFs
  red = "pca",             # dimensionality reduction to use
  group_vars = c("orig.ident")# , "Orientation")  # metadata for coloring
)


# Choose best PC/resolution combo from PDFs and then rerun on main object



# 4. Clustering

# FindNeighbors constructs a k-nearest neighbor (KNN) graph in the high-dimensional PCA space.
# It refines the edge weights based on shared local neighborhoods, producing a Shared Nearest Neighbor (SNN) graph.
Mouse_obj <- FindNeighbors(Mouse_obj, dims = 1:24)

# FindClusters is a separate modularity optimization step (typically Louvain) applied to the SNN graph.
# It groups cells into communities (clusters). This is an independent calculation from UMAP.
# - resolution: higher => more clusters
Mouse_obj <- FindClusters(Mouse_obj, resolution = 0.1)


DefaultAssay(Mouse_obj) <- "Xenium"

ImageDimPlot(Mouse_obj, fov = "X8fov", cols = Mouse_obj@misc$cluster_colors, size = 0.5, border.size = NA,
axes = TRUE, dark.background = TRUE) + DarkTheme()
ImageDimPlot(Mouse_obj, fov = "X8fov", molecules = "Slc17a7", nmols = 10000, alpha = 0.3, mols.cols = "red")



# 5. Non-linear Embedding: UMAP
# UMAP parameters:
# - dims: PCs to include (e.g., 1:24)
# Note: By default, RunUMAP independently computes its own nearest-neighbor graph in the high-dimensional 
# PCA space to estimate the manifold before mapping to the low-dimensional embedding. It does not use 
# the SNN graph created by FindNeighbors unless explicitly instructed via the nn.name parameter.

# Please note that there exist more than 1 implementation of the UMAP alogorithm.
# Seurat uses uwot (https://github.com/jlmelville/uwot) by default while there is a separate python implementation, 
# To run using umap.method="umap-learn", you must first install the umap-learn python package
# (e.g. via pip install umap-learn). For further details: https://github.com/lmcinnes/umap

# While the results are generally comparable there will be specific differences, as 
# shown here: https://github.com/satijalab/seurat/issues/2025. If you are trying to reproduce a 
# calculation of a lab member or something from the literature you may need to use the same implementation
# as the original calculation to ensure reproducibility.

# For an intuitive explanation of UMAP, see this StatQuest video: https://www.youtube.com/watch?v=eN0wFzBA4Sc
# paper on UMAP: https://arxiv.org/abs/1802.03426

# Here is a very useful playground with visualizations: https://pair-code.github.io/understanding-umap/

Mouse_obj <- RunUMAP(Mouse_obj, dims = 1:24)


# 6. Visualizations

# Standardize colour identifiers for each cluster to ensure consistency across graphs

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
DimPlot(Mouse_obj, #group.by = "seurat_clusters", 
        reduction = "umap",
        cols = Mouse_obj@misc$cluster_colors) + DarkTheme() + coord_fixed() + ggtitle("UMAP REDUCTION")


# You can plot your UMAP by metadata (ex. orig.ident, sex, stage, disease condition etc.)
#this is only after you have merged and/or integrated your data and added metadata columns
#DimPlot(Mouse_obj, reduction = "umap", cols = Mouse_obj@misc$cluster_colors) + DarkTheme() + coord_fixed() + ggtitle("UMAP REDUCTION grouped"



# You can Also cluster into independent objects to discern different processing techniques
#Mouse_Neighbors <- FindNeighbors(Mouse_obj, dims = 1:20)
#Mouse_Clustered <- FindClusters(Mouse_Neighbors, resolution = 0.1)


# ImageDimPlot(Mouse_Clustered, fov = "M1", cols = Mouse_obj@misc$cluster_colors, size = 0.5, border.size = NA,
#axes = TRUE, dark.background = TRUE) + DarkTheme()
# ImageDimPlot(Mouse_Clustered, molecules = "Slc17a7", nmols = 10000, alpha = 0.3, mols.cols = "red")

# Note: After running SCTransform, it sets the active assay to "SCT" (normalized data).
# For FeaturePlots() and VlnPlots() plot raw Xenium counts, temporarily switch the default assay.
# DefaultAssay(Mouse_Clustered) <- "Xenium"



# FeaturePlot: gene expression in UMAP space 
FeaturePlot(Mouse_obj, features = c("Aqp4"), label = TRUE,)
  
FeaturePlot(Mouse_obj, features = c("nCount_Xenium"), label = TRUE,)

VlnPlot(Mouse_obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size =0, log = TRUE)


# Gene expression in spatial coordinates
ImageFeaturePlot(Mouse_obj, fov = "X8fov", features = c("Aqp4", "Paqr5", "Trem2"))

#for single gene
ImageFeaturePlot(Mouse_obj, fov = "X8fov", features = ("Aqp4"))

ImageFeaturePlot(Mouse_obj, fov = "X8fov", features = ("nCount_Xenium"))

p <- ImageFeaturePlot(Mouse_obj, fov = "X8fov", features = "nCount_Xenium")

# Overwrite the default scale with a log10 transformation
p + scale_fill_gradient(low = "lightgrey", high = "firebrick", trans = "log10")

# Violin Plot of raw counts (optional: set log = TRUE for log scale)
VlnPlot(Mouse_obj, features = c("Slc17a7", "Gfap"), pt.size = 0, log = TRUE,)

# Summarize expression of selected genes across groups.
DotPlot(object = Mouse_obj, features = c("Aqp4", "Gfap", "Sla", "Paqr5", "Trem2"), dot.min  = 0.1,
        dot.scale= 6, group.by = "seurat_clusters") +
  scale_color_gradientn(colors = c("steelblue", "white", "firebrick")) +
  theme_classic() + ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),)



# Note: The above QC and Visualization is not just a linear pipeline, but should be
# treated as circular revision of filters and UMAP dims/neighbors based on the graphs 
# you see above.
# Eg: if a cluster is all just low-expression cells, its not very useful, so 
# you might want to increase nCount filter.
# Or maybe SCTransform does'nt function well with your dataset, so you use 
# log normalization etc. There could be many such scenarios to be mindful of.

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

