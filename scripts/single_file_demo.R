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

# Load in the RDS object
MAID_obj<-readRDS("M_70_Seurat_obj.rds")
# Update Seurat Objects to new structure for storing data/calculations. 
# For Seurat v3 objects, will validate object structure ensuring all keys 
# and feature names are formed properly.
MAID_obj <- UpdateSeuratObject(MAID_obj)


# Start by analyzing one sample at a time to understand the workflow and your data

# ------------------------------------------------------------------
# Single-Sample Analysis Workflow (QC, SCTransform, PCA, UMAP, Clustering)
# ------------------------------------------------------------------

## Preliminary analysis

# 1. Quality Control
# Visualize QC metrics

# Histogram 
# - nFeature_Xenium: number of unique transcript features detected per cell (i.e., how many genes/transcripts)
ggplot(MAID_obj@meta.data, aes(x = nFeature_Xenium)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  theme_classic() +
  scale_y_log10() +
  labs( #labels*
    title = "nFeature_Xenium",
    x = "Number of detected features",
    y = "Cell count"
  )

# - nCount_Xenium: total transcript count per cell (sum of all detected transcript counts)
# nCount_Xenium histogram
ggplot(MAID_obj@meta.data, aes(x = nCount_Xenium)) +
  geom_histogram(bins = 25, fill = "firebrick", color = "black") +
  theme_classic() +
  scale_y_log10() +
  labs(
    title = "nCount_Xenium",
    x = "Total transcript counts",
    y = "Cell count"
  )

#RidgePlot(MAID_obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2)



# Filter cells based on QC thresholds
#MAID_obj <- subset(MAID_obj, subset = nFeature_Xenium > 2 & nCount_Xenium > 0) 
# Note: to compare multiple QC filters for downstream analysis
# assign filtered results to new object names 
# (e.g., obj1filtered) instead of overwriting obj1, obj2 etc.
#MAID_obj_filtered <- subset(MAID_obj, subset = nFeature_Xenium > 0 & nCount_Xenium > 0)



#n_genes <- dim(MAID_obj)[1]
#n_cells <- dim(MAID_obj)[2]


# 2. Normalization and Feature Selection with SCTransform
# SCTransform performs normalization, variance stabilization, and identifies variable features
MAID_obj <- SCTransform(MAID_obj, assay = "Xenium")

#TBI <- SCTransform(TBI, assay = "Xenium")

#You can choose the normalization method and other details this way:
# 02_run_seurat.R lines 11-109
# Link to documentation ons SCTransform: https://satijalab.org/seurat/articles/sctransform_vignette.html
# Link to paper on SCTransform: https://link.springer.com/article/10.1186/s13059-021-02584-9

# function to plot residual variance vs gene expression 
residualVarPlot <- function(gene_var, xaxis = "gmean", max_resvar = 100, ntop = 25, annotate = F, pt_size = 1.1) {
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
  
  return(p)
}

# test normalization by plotting residual vs size of chart
# working on charting
gene_attr_maid <- SCTResults(MAID_obj, slot = "feature.attributes", assay = "SCT")
gene_attr_tbi <- SCTResults(MAID_obj, slot = "feature.attributes", assay = "SCT")

residualVarPlot(gene_attr_maid, max_resvar = 10, pt_size = 4)
residualVarPlot(gene_attr_tbi, max_resvar = 30, pt_size = 4)



# 3. Dimensionality Reduction: PCA
# PCA reduces high-dimensional data to principal components
# - features: use VariableFeatures from SCTransform
# - npcs: number of PCs to compute (e.g., 50)
# - Balance is Key: 
# - Including more PCs captures more variation in the data, but may introduce noise.
# - Including fewer PCs reduces noise but may risk missing important biological signals.

MAID_obj <- RunPCA(MAID_obj, features = VariableFeatures(MAID_obj), npcs = 50, verbose = FALSE)


# Option 1: Examine ElbowPlot to choose how many PCs to use downstream

# In this dataset, the elbow appears around PCs 12-15. 
# This suggests that most of the biologically relevant variation is captured within the first set of components.
# PC's after this point may be mostly noise. Lower PC dimensions may decrease this noise but could also miss
# biologically relevant information.
# PCA is about maximizing variance in the first components. However, 
# variance is a quantity that makes most sense for continuous and non-sparse data.
# if you dont see a clear elbow, you data might be unsuitable for an elbow plot examination
ElbowPlot(MAID_obj, ndims = 50)

# Option 2: PCA Heatmap
# - Use DimHeatmap() to visualize the top features (genes) driving each PC.
# - Adjust 'dims' to specify which PCs to examine; 'cells' controls how many cells to display per PC.
# - 'balanced = TRUE' scales positive and negative loadings equally.
DimHeatmap(MAID_obj, dims = 1:12, cells = 500, balanced = TRUE)



# 4. Non-linear Embedding: UMAP
# UMAP parameters:
# - dims: PCs to include (e.g., 1:20)
MAID_obj <- FindNeighbors(MAID_obj, dims = 1:24)
MAID_obj <- FindClusters(MAID_obj, resolution = 0.5)
MAID_obj <- RunUMAP(MAID_obj, dims = 1:24)

#MAID_UMAP_obj <- RunUMAP(MAID_obj, dims = 1:24)

#Plot UMAP (recommend adding PC and resol into the title to keep track of diff umaps)
DimPlot(MAID_obj, group.by = "seurat_clusters", 
        reduction = "umap",
        cols = "polychrome") + DarkTheme() + coord_fixed() + ggtitle("UMAP REDUCTION")


#Plot your UMAP by metadata (ex. orig.ident, sex, stage, disease condition etc.)
#this is only after you have merged and/or integrated your data and added metadata columns (Section 3 and 4)
#DimPlot(MAID_obj, reduction = "umap", cols = "polychrome") + DarkTheme() + coord_fixed() + ggtitle("UMAP REDUCTION grouped"



# 5. Clustering
# Identify clusters
# - resolution: higher => more clusters
MAID_Neighbors <- FindNeighbors(MAID_obj, dims = 1:20)
MAID_Clustered <- FindClusters(MAID_Neighbors, resolution = 0.1)


#ImageDimPlot(MAID_Clustered, fov = "M1", cols = "polychrome", size = 0.5, border.size = NA,
#axes = TRUE, dark.background = TRUE) + DarkTheme()
ImageDimPlot(MAID_Clustered, molecules = "SLC17A7", nmols = 10000, alpha = 0.3, mols.cols = "red")

# Note: After running SCTransform, it sets the active assay to "SCT" (normalized data).
# For FeaturePlots() and VlnPlots() plot raw Xenium counts, temporarily switch the default assay.
DefaultAssay(MAID_Clustered) <- "Xenium"
#FeaturePlot: gene expression in UMAP space 
FeaturePlot(MAID_obj, features = c("AQP4"), label = TRUE)

#Violin Plot of raw counts (optional: set log = TRUE for log scale)
VlnPlot(MAID_obj, features = c("SLC17A7", "GAD1", "GAD2"), pt.size = 0, log = TRUE)
#uses the polychrome colours
VlnPlot_scCustom(MAID_obj, features = c("SLC17A7", "GAD1", "GAD2"), pt.size = 0, log = TRUE)

DotPlot(object = MAID_obj, features = c("AQP4", "PAX6", "TGFB2"), dot.min  = 0.1,
        dot.scale= 6, group.by = "nFeature_Xenium") +
  scale_color_gradientn(colors = c("#0047AB", "white", "firebrick")) +
  theme_classic() + ggtitle("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))
ImageFeaturePlot(MAID_obj, features = c("AQP4", "PAX6", "TGFB2"))
#for single gene (works great)
ImageFeaturePlot(MAID_obj, fov = "M_70", features = ("SLC17A7"))





# ------------------------------------------------------------------
# Subsetting Data Based on Clusters of Interest
# ------------------------------------------------------------------
# If you are aiming to analyze a specific cell type or region you most likely need to subset your data
# Use a combination of marker gene expression and spatial location
# to identify the clusters you want to analyze further.
# you may need to perform several subsets to get to your cell type or region of interest based on how large your dataset is

###Several ways to subset (dependent on dataset)
##Option 1: Subset the object by clusters of interest (most common):
obj_sub <- subset(MAID_obj, idents = c(0, 2, 4, 6))

##Option 2: subset based on marker expression thresholds (e.g. cells that express a marker gene or combination of them)
# Example: keep cells where MarkerGene1 expression > 1 in the subset
#obj_sub <- subset(MAID_obj, subset = SLC17A7 > 1)

##Option 3: Crop your image down to a region of interest 
#Copy your original object
#obj1_crop <- MAID_obj

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
ImageDimPlot(obj_sub, fov = "M_70", cols = "polychrome", size = 0.5, border.size = NA,
             axes = TRUE, dark.background = TRUE) + DarkTheme()

##Once you have your subset using  any one of the options above, you can perform steps 1-7 as above
# Before re-running SCTransform on the refined subset, switch default assay
DefaultAssay(obj_sub) <- "Xenium"


#############SECTION 3: MERGING DATA##################################################

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
#merged_obj$Sex[merged_obj$orig.ident %in% c("F1"] <- "Female"

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
      plots <- list(DimPlot(dat_mod, shuffle = TRUE, label = TRUE) + coord_fixed())
      for (group_var in group_vars) {
        plots <- append(plots, list(DimPlot(dat_mod, group.by = group_var, shuffle = TRUE) + coord_fixed()))
      }
      combined_plot <- wrap_plots(plots, ncol = 1)
      plot_title <- paste("Dimensions:", max_dim, "| Resolution:", res)
      combined_plot <- combined_plot + plot_annotation(title = plot_title)
      file_name <- paste0("Dim_", max_dim, "_Res_", res, ".pdf")
      pdf_path <- file.path(output_dir, file_name)
      pdf(pdf_path, width = 8, height = 5 * length(plots))
      print(combined_plot)
      dev.off() # loop through clustering on each level of dimensions
    }
  }
}
# Run the automated test
dim_reso_test(
  dataset = merged_obj,
  dims_list = list(1:20, 1:30, 1:40, 1:50),  # test PC ranges
  resolutions = c(0.1, 0.4, 0.8),    # test resolutions
  # "Z:/Wellington Lab/Mehwish/Xenium_human_run1/B_HumanXenium_Run1_OG_rds"
  output_dir = "D:/work/Xenium/Output", # where to save PDFs
  red = "pca",             # dimensionality reduction to use
  group_vars = c("orig.ident", "Condition")  # metadata for coloring
)

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



