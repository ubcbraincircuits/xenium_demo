############SECTION 6: SINGLE-R ANALYSIS#######################

# Goal: Annotate Xenium (spatial) cells with cell-type labels using SingleR,
#       leveraging a single-cell RNA-seq reference.
#
# Overview of how it works:
# 1) Load a scRNA-seq reference (sc.ref) and one Xenium Seurat object (dat).
# 2) Normalize the scRNA-seq reference and (optionally) pseudobulk it by a
#    label (here: "majorcelltype") to create a robust, denoised reference matrix.
# 3) For each sample in the Xenium object (dat$orig.ident), normalize the Xenium
#    assay, intersect genes with the reference, and run SingleR to predict labels.
# 4) Save per-sample SingleR results and write the predicted labels back to the
#    full objects’ metadata (both raw and "pruned" labels for stricter calls).
# by larissa June, 2025
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
sc <- readRDS("yourscRNAseqobj.rds")
# Xenium data 
ST <- readRDS("merged_obj.rds")


# assign scRNAseq as your reference 
sc.ref <- sc

# make dat your xenium obj 
dat <- cortical_ST


# quick checks on reference annotations and assays
unique(sc.ref$majorcelltype)
Assays(sc.ref)

# Set the output directory for saving per-sample SingleR result RDS files
output_dir <- 'X:/Cembrowski Lab/pathtoyouroutputfolder/'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

## Prepare the reference expression matrix
# Log-Normalize sc.ref (SingleR expects normalized/log-scale expression)
sc.ref <- NormalizeData(sc.ref, assay = "RNA", verbose = FALSE)

# Pseudobulk the single-cell reference (recommended: faster and less noisy than per-cell)
pseudobulk <- AggregateExpression(
  sc.ref,
  group.by = "majorcelltype",  # labels defining reference classes
  assays = "RNA",
  slot = "data", 
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
  query_expr <- GetAssayData(sample_obj, slot = "data", assay = "Xenium")
  
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

#visualize using methods discussed above
#ensure to 'group.by = SingleR_pruned_midfine or SingleR_label_midfine
