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


# Add SingleR labels
singler_results <- readRDS("Output/SingleR_ST_results.rds")
Mouse_obj$CellSubtype <- singler_results$SingleR_pruned_midfine

Mouse_obj$condition <- NA

Mouse_obj$condition[Mouse_obj$orig.ident %in% c("X2_fov", "X3_fov")] <- "A"
Mouse_obj$condition[Mouse_obj$orig.ident %in% c("X4_fov", "X5_fov")] <- "B"

# Keep only cells with assigned condition and non-missing cell subtype
Mouse_obj <- subset(
  Mouse_obj,
  subset = !is.na(condition) & !is.na(CellSubtype)
)

# Check metadata
table(Mouse_obj$orig.ident, Mouse_obj$condition)
table(Mouse_obj$CellSubtype, Mouse_obj$condition)

# Add SingleR labels
singler_results <- readRDS("Output/SingleR_ST_results.rds")
Mouse_obj$CellSubtype <- singler_results$SingleR_pruned_midfine

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
pseudo <- AggregateExpression(
  Mouse_obj,
  assays = "Xenium",
  return.seurat = TRUE,
  group.by = c("condition", "orig.ident", "CellSubtype")
)

# Create identity for same-celltype comparison across conditions
pseudo$celltype.condition <- paste(pseudo$CellSubtype, pseudo$condition, sep = "_")
Idents(pseudo) <- "celltype.condition"

# Inspect available groups
table(Idents(pseudo))
head(pseudo@meta.data)

# Note, generally DESeq2 and FindMarkers needs a minimum of 3 samples per group
# Below is just a demonstration of how to do pseudobulk DEG in xenium data, but
# it is highly recommended to have more samples for actual analysis

# Run DE between pseudobulk groups
pb_markers_astro <- FindMarkers(
  object = pseudo,
  ident.1 = "Astro_A",
  ident.2 = "Astro_B",
  assay = "Xenium",
  test.use = "DESeq2",
  min.cells.group = 2
)

head(pb_markers_astro)

pb_markers_l4 <- FindMarkers(
  object = pseudo,
  ident.1 = "L4_A",
  ident.2 = "L4_B",
  assay = "Xenium",
  test.use = "DESeq2",
  min.cells.group = 2
)


 ###Visualize marker genes using appropriate plots for pseudobulked data

##Boxplots 
# Boxplots are preferred for pseudobulked data since each point 
# represents an aggregated group (e.g., sample),
# and violin plots can misrepresent distribution when cell-level resolution is lost.

Idents(pseudo) <- "celltype.condition"
DefaultAssay(pseudo) <- "Xenium"

plot_pb_box <- function(pseudo, gene, celltype = "Astro", group1 = "A", group2 = "B") {
  ids <- c(paste(celltype, group1, sep = "_"),
           paste(celltype, group2, sep = "_"))
  
  df <- FetchData(
    object = pseudo,
    vars = c("celltype.condition", "condition", "CellSubtype", gene)
  )
  
  colnames(df)[ncol(df)] <- "expr"
  df <- df %>% filter(celltype.condition %in% ids)
  
  ggplot(df, aes(x = celltype.condition, y = expr, color = condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.45, linewidth = 0.35) +
    geom_point(size = 2.8, position = position_jitter(width = 0.08, height = 0)) +
    theme_classic() +
    labs(
      title = paste0(gene, " pseudobulk expression in ", celltype),
      x = NULL,
      y = "Aggregated expression"
    ) +
    scale_color_manual(values = c("A" = "#E64B35", "B" = "#4DBBD5"))
}

# Example
p_box_aqp4 <- plot_pb_box(pseudo, gene = "Aqp4", celltype = "Astro")
p_box_aqp4

genes_to_plot <- c("Aqp4", "Paqr5", "Trem2")

boxplots <- lapply(genes_to_plot, function(g) {
  plot_pb_box(pseudo, gene = g, celltype = "L4")
})

boxplots[[1]]
boxplots[[2]]
boxplots[[3]]



Idents(pseudo) <- "Celltype" #Switch to correct ident


# note that there are'nt enough samples to make a representative plot here
VlnPlot(
  object = pseudo,
  features = "Aqp4",
  assay = "Xenium",
  idents = c("Astro_A", "Astro_B"),
  pt.size = 1,
  adjust = 0
) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    linewidth = 0.3,
    fill = NA
  ) +
  theme_classic()


plot_pb_volcano <- function(markers_df,
                            title = "Pseudobulk DE",
                            x_thresh = 1,
                            p_thresh = 0.05,
                            label_top = 10,
                            up_label = "Up in A",
                            down_label = "Up in B") {
  
  df <- markers_df %>%
    rownames_to_column("gene") %>%
    mutate(
      p_val_adj = ifelse(is.na(p_val_adj), 1, p_val_adj),
      p_plot = pmax(p_val_adj, 1e-300),
      status = case_when(
        p_val_adj < p_thresh & avg_log2FC >= x_thresh  ~ up_label,
        p_val_adj < p_thresh & avg_log2FC <= -x_thresh ~ down_label,
        TRUE ~ "NS"
      )
    )
  
  lab_df <- df %>%
    filter(status != "NS") %>%
    arrange(p_val_adj) %>%
    slice_head(n = label_top)
  
  ggplot(df, aes(x = avg_log2FC, y = -log10(p_plot), color = status)) +
    geom_point(alpha = 0.8, size = 1.8) +
    geom_vline(xintercept = c(-x_thresh, x_thresh), linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", linewidth = 0.3) +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = gene),
      size = 3,
      max.overlaps = 20
    ) +
    scale_color_manual(values = c(
      "NS" = "grey70",
      up_label = "#E64B35",
      down_label = "#4DBBD5"
    )) +
    theme_classic() +
    labs(
      title = title,
      x = "avg_log2FC",
      y = "-log10(adjusted p-value)",
      color = NULL
    )
}

# Example with your DE result
p_volcano_astro <- plot_pb_volcano(
  pb_markers_astro,
  title = "Astro pseudobulk: A vs B",
  x_thresh = 1,
  p_thresh = 0.05,
  label_top = 12,
  up_label = "Higher in A",
  down_label = "Higher in B"
)

p_volcano_astro
