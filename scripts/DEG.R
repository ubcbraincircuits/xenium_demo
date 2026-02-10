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

MAID_obj = SCTransform(MAID_obj, assay = "Xenium")


MAID_obj <- RunPCA(MAID_obj, features = VariableFeatures(MAID_obj), npcs = 50, verbose = FALSE)
MAID_obj <- FindNeighbors(MAID_obj, dims = 1:20)
MAID_obj <- FindClusters(MAID_obj, resolution = 0.1)

# Assign global color pallete
global_clusters <- levels(Idents(MAID_obj))
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
MAID_obj@misc$cluster_colors <- cluster_colors

# Ensure Idents is a factor with consistent levels
Idents(MAID_obj) <- factor(Idents(MAID_obj), levels = global_clusters)

#############SECTION 5: DIFFERENTIAL GENE EXPRESSION (DEG) ANALYSIS#############################

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
######### 1. SINGLE-CELL LEVEL DEG USING SCTransform #########

# Function to prep SCT assay for differential expression testing
prep_FindMarkers <- function(obj, num_slices = length(obj@assays$SCT@SCTModel.list)) {
  # Set the raw UMI assay (used during SCTransform modeling) for each slice
  for (i in 1:num_slices) {
    slot(object = obj@assays$SCT@SCTModel.list[[i]], name = "umi.assay") <- "Xenium"
  }
  # Prepares the SCT assay for DE testing by computing required offsets
  obj <- PrepSCTFindMarkers(obj)
  return(obj)
}

# Apply the prep function to your Seurat object (using MAID_obj as example)
MAID_obj <- prep_FindMarkers(MAID_obj, num_slices = 1)

# Run DE across all clusters using SCT-normalized expression
sc_markers <- FindAllMarkers(
  object          = MAID_obj,
  assay           = "SCT",           # Use SCT-normalized values
  only.pos        = TRUE) #Change to FALSE for both 'upregulated' and 'downregulated' genes


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


DefaultAssay(MAID_obj) <- "SCT"


VlnPlot(
  MAID_obj,
  features = c("nFeature_Xenium", "nCount_Xenium"),
  group.by = "seurat_clusters",
  ncol = 2,
  pt.size = 0,
  cols = MAID_obj@misc$cluster_colors
) +
  theme_classic()

VlnPlot(
  MAID_obj,
  features = "AQP4",
  group.by = "seurat_clusters",
  pt.size = 0,
  cols = MAID_obj@misc$cluster_colors
) +
  theme_classic()


FeaturePlot(
  MAID_obj,
  features = c("SLC17A6", "GAD1"),
  #assay = "SCT",
  ncol = 3,
  min.cutoff = "q10",
  max.cutoff = "q90",
)



######### 2. PSEUDOBULK LEVEL DEG ANALYSIS ##########################################
# Reloading MAID_obj

# Aggregate raw counts by sample and cell subtype (or other metadata grouping)
pseudo <- AggregateExpression(MAID_obj, assays = "Xenium", return.seurat = TRUE,
                              group.by      = c("orig.ident", "CellSubtype"))

# Optional: check identities of resulting pseudobulk groups
table(Idents(pseudo))  # Should show sample × cell subtype combinations

# Run DE between pseudobulk groups
pb_markers <- FindMarkers(object = pseudo, ident.1 = "CellType1",      
                          ident.2 = "CellType2", assay = "Xenium", test.use = "DESeq2")


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