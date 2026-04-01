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

# First, once again lets clear old objects
rm(list = ls())
gc()


#######SECTION: LOAD IN RDS OBJECTS, PERFORM QUALITY CONTROL, NORMALIZE DATA, DIMENTIONALITY REDUCTION AND CLUSTERING###########
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

Mouse_obj = SCTransform(Mouse_obj, assay = "Xenium")


Mouse_obj <- RunPCA(Mouse_obj, features = VariableFeatures(Mouse_obj), npcs = 50, verbose = FALSE)
Mouse_obj <- FindNeighbors(Mouse_obj, dims = 1:24)
Mouse_obj <- FindClusters(Mouse_obj, resolution = 0.1)

Mouse_obj <- RunUMAP(Mouse_obj, dims = 1:24)

# Assign global color pallete
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

#############SECTION: DIFFERENTIAL GENE EXPRESSION (DEG) ANALYSIS #############################

# There are two primary approaches to performing DEG analysis, depending on the 
# resolution and robustness you require:
#
# 1. Single-cell DE (FindMarkers / FindAllMarkers):
#    Best for detecting gene differences between clusters or conditions at 
#    single-cell resolution, assuming you have a large enough number of cells 
#    and minimal technical noise. This approach uses normalized data (e.g., SCT).
#
# 2. Pseudobulk DE:
#    Preferred for group-level comparisons (e.g., per sample, region, or 
#    condition) when you want to reduce cell-level variability or use robust 
#    statistical frameworks like DESeq2.

###############################################################################################


######### 1. SINGLE-CELL LEVEL DEG USING SCTransform #########

# For merged objects with multiple SCT models, Seurat recommends running
# PrepSCTFindMarkers() before DE on the SCT assay.
# Documentation:
# https://satijalab.org/seurat/reference/prepsctfindmarkers
#
# This function standardizes the SCT representation so marker testing is more
# comparable across merged pieces / slices. 
prep_FindMarkers <- function(obj, num_slices = length(obj@assays$SCT@SCTModel.list)) {
  # Set the raw UMI assay (used during SCTransform modeling) for each slice
  for (i in 1:num_slices) {
    slot(object = obj@assays$SCT@SCTModel.list[[i]], name = "umi.assay") <- "Xenium"
  }
  # Prepares the SCT assay for DE testing by computing required offsets
  obj <- PrepSCTFindMarkers(obj)
  return(obj)
}

# Apply the prep function to your Seurat object (using Mouse_obj as example)
# For this single-object example, num_slices = 1 is sufficient.
Mouse_obj <- prep_FindMarkers(Mouse_obj, num_slices = 1)

# Run DE across all clusters using SCT-normalized expression
# FindAllMarkers identifies marker genes for each identity class (cluster)
# relative to all other cells.
# Documentation: https://satijalab.org/seurat/reference/findallmarkers
sc_markers <- FindAllMarkers(
  object          = Mouse_obj,
  assay           = "SCT",# Use SCT-normalized values
  only.pos        = TRUE) # Change to FALSE for both 'upregulated' and 'downregulated' genes


# After identifying (DEGs) with FindAllMarkers(),
# it is useful to apply additional filters to retain genes that are not
# only statistically significant, but also biologically meaningful.
#
# This is an example code that filters the DEGs:
# - Are expressed in at least 30% of cells in the cluster of interest (pct.1 > 0.3),
# - And either:
#   a) Are expressed in fewer than 30% of cells in other clusters (pct.2 < 0.3),
#   b) OR show a strong fold change (abs(log2FC) > 1).
#
# This results in a more focused and robust set of cluster-enriched marker genes.

# Filter DEGs 
sc_markers <- subset(sc_markers, pct.1 > 0.3 & (pct.2 < 0.3 | abs(avg_log2FC) > 1))

#Visualize marker genes using using FeaturePlots, VlnPlots, DotPlots as described above. 


# Extract the single best marker per cluster based on average log2 fold change
# to quickly visualize and verify
top1_markers <- sc_markers %>%
  group_by(cluster) %>%
  slice_max(n = 1, order_by = avg_log2FC)

# Create a clean character vector of these top genes to feed into plotting functions
top_genes <- unique(top1_markers$gene)

DefaultAssay(Mouse_obj) <- "SCT"

# percentage of cells expressing a gene and the average expression level across clusters
DotPlot(
  Mouse_obj, 
  features = top_genes, 
  group.by = "seurat_clusters"
) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Show the expression distribution of a gene across clusters
VlnPlot(
  Mouse_obj,
  features = top_genes[1:3], # Plucking the first 3 top markers 
  group.by = "seurat_clusters",
  pt.size = 0,
  cols = Mouse_obj@misc$cluster_colors
) +
  theme_classic()

VlnPlot(
  Mouse_obj,
  features = "Aqp4",
  group.by = "seurat_clusters",
  pt.size = 0,
  cols = Mouse_obj@misc$cluster_colors
) +
  theme_classic()

FeaturePlot(
  Mouse_obj,
  features = top_genes[1:3], 
  ncol = 3,
  min.cutoff = "q10",
  max.cutoff = "q90",
  reduction = "umap",
)

DimPlot(
  Mouse_obj,
  reduction = "umap",
  label = TRUE,        # Overlays the cluster names/numbers on the plot
  repel = TRUE,        # Prevents text labels from overlapping
  pt.size = 0.5
) +
  ggtitle("Clusters") +
  theme_classic()


FeaturePlot(
  Mouse_obj,
  features = top_genes[1],
  min.cutoff = "q10",
  max.cutoff = "q90",
  reduction = "umap",
)


ImageFeaturePlot(Mouse_obj, fov = "X8fov", features = c("Aqp4", "Paqr5", "Trem2"))



######### 2. PSEUDOBULK LEVEL DEG ANALYSIS ##########################################
# First, once again lets clear old objects
rm(list = ls())
gc()

# let's reload mouse object and apply prepossessing steps to it
#Mouse_obj<-readRDS("Region_2_right-obj.rds")
Mouse_obj<-readRDS("merged_raw_metadata_obj.rds")
# Update Seurat Objects to new structure for storing data/calculations. 
# For Seurat v3 objects, will validate object structure ensuring all keys 
# and feature names are formed properly.
Mouse_obj <- UpdateSeuratObject(Mouse_obj)

Mouse_obj = SCTransform(Mouse_obj, assay = "Xenium")


Mouse_obj <- RunPCA(Mouse_obj, features = VariableFeatures(Mouse_obj), npcs = 50, verbose = FALSE)
Mouse_obj <- FindNeighbors(Mouse_obj, dims = 1:24)
Mouse_obj <- FindClusters(Mouse_obj, resolution = 0.1)

Mouse_obj <- RunUMAP(Mouse_obj, dims = 1:24)

# get saved SingleR results from previous script
singler_results <-readRDS("Output/SingleR_ST_results.rds")
Mouse_obj$CellSubtype <- singler_results$SingleR_pruned_midfine
#Mouse_obj$CellSubtype <- singler_results$pruned.labels

# Assign global color pallete
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


#### Start of Psuedobulk specific code ###


# Aggregate raw counts by sample and cell subtype (or other metadata grouping)
pseudo <- AggregateExpression(Mouse_obj, assays = "Xenium", return.seurat = TRUE,
                              group.by      = c("orig.ident")) #, "CellSubtype"))

# Extract the cell type from the pseudobulk column names 
# (Assuming the format is orig.ident_CellSubtype, so we grab the second part)
#pseudo$CellType <- sapply(strsplit(colnames(pseudo), "_"), `[`, 2)

# Optional: check identities of resulting pseudobulk groups
table(Idents(pseudo))  # Should show sample × cell subtype combinations

# Run DE between pseudobulk groups
pb_markers <- FindMarkers(object = pseudo, ident.1 = "X2-fov_Astro",      
                          ident.2 = "X2-fov_Macrophage", assay = "Xenium", test.use = "wilcox_limma") #test.use = "DESeq2")


 ###Visualize marker genes using appropriate plots for pseudobulked data

##Boxplots 
# Boxplots are preferred for pseudobulked data since each point 
# represents an aggregated group (e.g., sample),
# and violin plots can misrepresent distribution when cell-level resolution is lost.

Idents(pseudo) <- "Celltype" #Switch to correct ident

VlnPlot(object   = pseudo, features = "Gene", assay = "Xenium",
        idents   = c("CellType1", "CellType2"), pt.size  = 0,
        adjust   = 0) + 
  geom_boxplot(
    outlier.size   = 0,
    outlier.stroke = 0,
    lwd = 0.3)


# Volcano plot showing differential expression between two pseudobulked groups.
VolcanoPlot(pseudo, 
            ident.1 = "CellType1",
            ident.2 = "CellType2", 
            y.threshold = 1, x.threshold = 1, log.base = "2") + #adjust thresholds as needed
  theme_classic()
