############SECTION 6: SINGLE-R ANALYSIS#######################

# Goal: Annotate Xenium (spatial) cells with cell-type labels using SingleR,
#       leveraging a single-cell RNA-seq reference.
#
# Overview of how it works:
# 1) Load a scRNA-seq reference (sc.ref) and one Xenium Seurat object (dat).
# 2) Normalize the scRNA-seq reference and (optionally) pseudobulk it by a
#    label (here: "subclass") to create a robust, denoised reference matrix.
# 3) For each sample in the Xenium object (dat$orig.ident), normalize the Xenium
#    assay, intersect genes with the reference, and run SingleR to predict labels.
# 4) Save per-sample SingleR results and write the predicted labels back to the
#    full objects’ metadata (both raw and "pruned" labels for stricter calls).
# by larissa June, 2025


# We choose to use SingleR over other cell type annotation methods as it is fast
# and accurate in comparison to alternatives as laid out in this paper:
# https://link.springer.com/article/10.1186/s12859-025-06044-0


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



# Set the working directory to the folder containing your saved RDS objects
setwd("D:/work/Github_demo/xenium_demo/Data")

# read in single-cell and xenium data 
# single cell data (reference)
#sc <- readRDS("allen_cortex.rds")
# Xenium data 
ST <- readRDS("merged_obj.rds")

# Normalize and run reductions   Not needed because normalization done in loop!!
#ST <- SCTransform(ST, assay = "Xenium")

#ST <- RunPCA(ST, features = VariableFeatures(ST), npcs = 50, verbose = FALSE)

#ST <- FindNeighbors(ST, dims = 1:30)
#ST <- FindClusters(ST, resolution = 0.5)
#ST <- RunUMAP(ST, dims = 1:30)


# assign scRNAseq as your reference 
#sc.ref <- sc

# make dat your xenium obj 
dat <- ST #cortical_ST

# single cell data (reference)
sc.ref <- readRDS("allen_cortex.rds") 
# from this link: https://satijalab.org/seurat/archive/v3.2/spatial_vignette.html 
# and this data: https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1


# Xenium data
#dat <- readRDS("merged_obj.rds") 

# Run garbage collection to free up temporary RAM overhead
gc()


# quick checks on reference annotations and assays
#unique(sc.ref$majorcelltype) # Allen_cortex reference uses subclass
unique(sc.ref$subclass)
Assays(sc.ref)

# Set the output directory for saving per-sample SingleR result RDS files
#output_dir <- 'X:/Cembrowski Lab/pathtoyouroutputfolder/'
output_dir <- 'D:/work/Github_demo/xenium_demo/Data/Output'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Prepare the reference expression matrix
# Log-Normalize sc.ref (SingleR expects normalized/log-scale expression)
sc.ref <- NormalizeData(sc.ref, assay = "RNA", verbose = FALSE)

# Pseudobulk the single-cell reference (recommended: faster and less noisy than per-cell)
pseudobulk <- AggregateExpression(
  sc.ref,
  #group.by = "majorcelltype",  # labels defining reference classes
  group.by = "subclass",
  assays = "RNA",
  #slot = "data",
  #layer = "data", 
  return.seurat = FALSE)

# Matrix of genes x cell types (normalized counts for SingleR "ref")
ref_expr <- pseudobulk$RNA

# Identify samples in the Xenium dataset; loop runs per sample for reproducibility
sample_names <- unique(dat$orig.ident)

# Pre-allocate holders for all-cell predictions (labels and pruned.labels)
all_labels <- rep(NA, ncol(dat))
all_pruned <- rep(NA, ncol(dat))
names(all_labels) <- colnames(dat)
names(all_pruned) <- colnames(dat)

# run loop to execute SingleR on each sample in dat
for (sample_name in sample_names) {
  message(paste0("Running SingleR for sample: ", sample_name))
  
  # Subset Seurat object for current sample
  sample_obj <- subset(dat, subset = orig.ident == sample_name)
  
  # Normalize Xenium assay (ensure query is on log-normalized scale)
  sample_obj <- NormalizeData(sample_obj, assay = "Xenium", verbose = FALSE)
  
  # Extract log-normalized expression for the query (Xenium)
  query_expr <- GetAssayData(sample_obj, layer = "data", assay = "Xenium") # slot = "data"
  
  # Intersect genes between reference and query to match features
  # (Double subset avoids control probes like blank_codewords)
  common_genes <- intersect(rownames(ref_expr), rownames(query_expr))
  ref_sub <- ref_expr[common_genes, , drop = FALSE]
  query_sub <- query_expr[common_genes, , drop = FALSE]
  
  # Run SingleR using pseudobulked reference; labels are the ref column names
  singleR_results <- SingleR(
    test = as.matrix(query_sub),
    ref = as.matrix(ref_sub),
    labels = colnames(ref_expr)
  )
  
  # Save per-sample SingleR results for audit/reuse
  saveRDS(singleR_results, file = file.path(output_dir, paste0("SingleR_", sample_name, "_results.rds")))
  
  # Map predictions back to the full 'dat' metadata using barcodes
  # pruned.labels = higher-confidence subset of calls
  sample_cells <- colnames(sample_obj)
  all_labels[sample_cells] <- singleR_results$labels
  all_pruned[sample_cells] <- singleR_results$pruned.labels
}

# Add to metadata of the full dat object(s) 

ST$SingleR_label_midfine <- all_labels
ST$SingleR_pruned_midfine <- all_pruned

#labeled_merged_obj <- AddMetaData(dat, metadata = all_labels, col.name = "CellSubtype_SingleR")

saveRDS(ST, file = file.path(output_dir, paste0("SingleR_ST_results.rds")))


#visualize using methods discussed above
#ensure to 'group.by = SingleR_pruned_midfine or SingleR_label_midfine

# Run garbage collection to free up temporary RAM overhead
gc()

## Note adjust fov to whatever you want to visualize!!
# find fov's using:

Images(ST) # or use names(obj@images)


# Visualize the high-confidence annotations on one specific FOV
ImageDimPlot(ST, 
             fov = "X2fov",  # Change this to whichever FOV you want to inspect
             group.by = "SingleR_pruned_midfine", 
             size = 0.5, 
             dark.background = TRUE)

# Create a new column just to highlight unannotated cells
ST$is_unannotated <- is.na(ST$SingleR_pruned_midfine)
ImageDimPlot(ST, 
             fov = "X2fov", 
             group.by = "is_unannotated", 
             cols = c("FALSE" = "grey", "TRUE" = "red"),
             size = 0.5, 
             dark.background = TRUE)

# Create a proportional bar chart comparing left vs right, or region vs region
ggplot(ST@meta.data, aes(x = orig.ident, fill = SingleR_pruned_midfine)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(x = "Cortex Region / FOV", y = "Proportion of Cells", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Normalize, cluster and Dim reduction
ST <- FindVariableFeatures(ST, selection.method = "vst", nfeatures = 2000)

