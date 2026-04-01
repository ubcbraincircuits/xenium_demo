############ SECTION: SINGLE-R ANALYSIS (Single FOV) #######################

# Goal: Annotate Xenium (spatial) cells with cell-type labels using SingleR,
# leveraging a single-cell RNA-seq reference.
#
# Overview of how it works:
# 1) Load a scRNA-seq reference (sc.ref) and one Xenium Seurat object (Mouse_obj).
# 2) Normalize the scRNA-seq reference and (optionally) pseudobulk it by a
#    label (here: "subclass") to create a robust, denoised reference matrix.
# 3) Normalize the Xenium assay, intersect genes with the reference, and run 
#    SingleR to predict labels on a single FOV object.
# 4) Save SingleR results and write the predicted labels back to the
#    objects’ metadata (both raw and "pruned" labels for stricter calls).
# Modified code written by larissa June, 2025

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

# Set the working directory to the folder containing your saved RDS objects
setwd("D:/work/Github_demo/xenium_demo/Data")

# First, we get a reference scRNA-seq dataset. Here a dataset of ~14,000 adult mouse 
# cortical cell taxonomy from the Allen Institute, generated with the SMART-Seq2 protocol.
# link: https://www.nature.com/articles/nn.4216
# You should have downloaded this alongside the other RDS files from OSF.
sc.ref <- readRDS("allen_cortex.rds") 
# Can also be accessed here: https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1

# Xenium data for a single FOV/sample. Our Query data
Mouse_obj <- readRDS("Region_2_right-obj.rds")
Mouse_obj <- UpdateSeuratObject(Mouse_obj)

# Run garbage collection to free up temporary RAM overhead
gc()

# Set the output directory for saving SingleR result RDS files
output_dir <- 'D:/work/Github_demo/xenium_demo/Data/Output'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Prepare the reference expression matrix

# Log-Normalize sc.ref (SingleR expects normalized/log-scale expression)
# Note we don't use SCTransform here, log-normalization is also an accessible normalization
# function for datasets, where SCTransform is not neccessary/applicable
sc.ref <- NormalizeData(sc.ref, assay = "RNA", verbose = FALSE)

# Pseudobulk the single-cell reference.
# It combines the expression profiles of all cells belonging to a given 
# cell type (here, defined by the "subclass" metadata column).
# It creates a cleaner, denoised average reference profile, minimizing single-cell dropouts.
# And it significantly speeds up SingleR computation since it compares against a handful of 
# cell types rather than thousands of individual reference cells.
pseudobulk <- AggregateExpression(
  sc.ref,
  group.by = "subclass",
  assays = "RNA",
  return.seurat = FALSE
)

# Matrix of genes x cell types (normalized counts for SingleR "ref")
# pseudobulk$RNA returns the aggregated matrix
ref_expr <- pseudobulk$RNA

# Normalize Xenium assay (ensure query is on log-normalized scale, just like ref)
Mouse_obj <- NormalizeData(Mouse_obj, assay = "Xenium", verbose = FALSE)

# Extract log-normalized expression for the query
# GetAssayData retrieves the specified data layer from a Seurat object.
query_expr <- GetAssayData(Mouse_obj, layer = "data", assay = "Xenium")

# Intersect genes between reference and query to match features
# A classifier can only evaluate features (genes) that exist in both datasets
# (Double subset avoids control probes like blank_codewords)
common_genes <- intersect(rownames(ref_expr), rownames(query_expr))

# Subset both matrices to only the overlapping features
ref_sub <- ref_expr[common_genes, , drop = FALSE]
query_sub <- query_expr[common_genes, , drop = FALSE]

message("Running SingleR for query...")

# Run SingleR using pseudobulked reference
# SingleR is a label transfer technique that calculates the Spearman correlation
# between each single cell in the query and the reference profiles. Then It fine-tunes
# these predictions by comparing only the most informative marker genes between
# closely correlated cell types.

# More detailed explanation can be found here: https://bioconductor.org/books/release/SingleRBook/introduction.html

# It is held as a fast and accurate label-transfer technique in comparison to 
# alternatives accessible to seurat. 
# As explained in this paper: https://link.springer.com/article/10.1186/s12859-025-06044-0
# labels here are the corresponding cell type names for the columns of the ref matrix
singleR_results <- SingleR(
  test = as.matrix(query_sub),
  ref = as.matrix(ref_sub),
  labels = colnames(ref_expr)
)

# Save SingleR results for audit/reuse
saveRDS(singleR_results, file = file.path(output_dir, "SOLO_SingleR_X2_FOV_results.rds"))

# Map predictions back to the full 'Mouse_obj' metadata using barcodes
# labels are the raw, best-guess prediction for every cell.
# pruned.labels are a higher-confidence subset of calls. Cells that have ambiguous or 
# low-confidence correlations are set to NA (pruned), making your final annotations more robust
Mouse_obj$SingleR_label_midfine <- singleR_results$labels
Mouse_obj$SingleR_pruned_midfine <- singleR_results$pruned.labels

# Optional: Save the updated Seurat object
# saveRDS(Mouse_obj, file = file.path(output_dir, "SingleR_Mouse_obj_results.rds"))

# Run garbage collection to free up temporary RAM overhead
gc()

## Visualization
# Find the name of the FOV dynamically
fov_name <- Images(Mouse_obj)[1] 

# Visualize the high-confidence annotations on the specific FOV
ImageDimPlot(Mouse_obj, 
             fov = fov_name, 
             group.by = "SingleR_pruned_midfine", 
             size = 0.5, 
             dark.background = TRUE)

# Create a new column just to highlight unannotated cells in pruned labels
Mouse_obj$is_unannotated <- is.na(Mouse_obj$SingleR_pruned_midfine)
ImageDimPlot(Mouse_obj, 
             fov = fov_name, 
             group.by = "is_unannotated", 
             cols = c("FALSE" = "grey", "TRUE" = "red"),
             size = 0.5, ++
             dark.background = TRUE)

# Create a proportional bar chart of cell types for the single FOV
ggplot(Mouse_obj@meta.data, aes(x = orig.ident, fill = SingleR_pruned_midfine)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(x = "FOV", y = "Proportion of Cells", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Normalize, cluster and Dim reduction ...
Mouse_obj <- FindVariableFeatures(Mouse_obj, selection.method = "vst", nfeatures = 2000)