# Xenium Spatial Transcriptomics Analysis Scripts Workflow

Please go through the scripts in the numerical order provided below. Each script builds upon concepts and files generated in the previous steps.

## Step 0: Setting Up Your Environment
**File:** `0_Load_Packages.R`

Before analyzing any data, you must ensure your R environment is correctly configured. 
- **What it does:** This script safely installs and loads all the necessary packages for the entire workflow (e.g., `Seurat`, `ggplot2`, `SingleR`, `DESeq2`). It uses `BiocManager` to fetch tools from both CRAN and Bioconductor.

## Step 1: Signle Sample Basics
**File:** `1_Single_File_Demo.R`

We start analysis by looking at just one Field of View (FOV) / sample. This makes it easier to understand the data structure before dealing with multi-sample datasets.
- **What it does:** Walks you through loading a `.rds` Seurat object and exploring its metadata and spatial coordinates. 
- **Key Concepts Covered:** 
  - **Quality Control (QC):** Filtering out poor-quality cells based on transcript counts.
  - **Normalization:** Using `SCTransform` to remove technical noise while preserving true biological variance.
  - **Dimensionality Reduction & Clustering:** Running PCA, selecting dimensions via Elbow plots, grouping cells into clusters, and visualizing them with UMAP.
  - **Visualization:** Plotting gene expression directly onto the spatial tissue map and/or UMAP space

## Step 2: Combining Multiple Samples
**File:** `2_Merging_Demo.R`

Here, we combine multiple spatial samples into one dataset for convenience in downstream analysis.
- **What it does:** Demonstrates how to merge multiple FOVs into a single Seurat object using simple Merge and/or Harmony integration.
- **Key Concepts Covered:**
  - **Merging:** Stitching together raw `.rds` objects.
  - **Batch Correction (Harmony):** active "Integration" using Harmony.

## Step 3: Cell Type Label Transfer
**File:** `3_SingleR.R`

Here, we give cells a cell-type annotation (e.g., Astrocytes, Neurons, Microglia).
- **What it does:** Uses `SingleR` to automatically predict cell types in your spatial dataset by comparing it to a scRNA-seq reference dataset.
- **Key Concepts Covered:**
  - **Reference Pseudobulking:** Aggregating the scRNA-seq reference to create a clean, denoised profile for faster and more accurate matching.
  - **Label Transfer:** Mapping the reference labels onto both a single FOV and your fully merged dataset.
  - **Confidence Pruning:** Dropping ambiguous cells to keep only high-confidence annotations (`SingleR_pruned_midfine`).

## Step 4: DEG
**File:** `4_Differential_expression.R`

With your cells clustered and identified, the final step is to find out which genes are significantly upregulated or downregulated between different groups.
- **What it does:** Different statistical approaches to finding Differentially Expressed Genes (DEGs).
- **Key Concepts Covered:**
  - **Single-Cell DEGs:** Using `FindAllMarkers` on `SCTransform` normalized data to find genes that define specific clusters.
  - **Pseudobulk DEGs:** We use `AggregateExpression` to sum counts by replicate, and test them using `DESeq2`.

> NOTE: This walkthrough is meant as a demonstration of some of the things you can do with Seurat and Xenium data. 
> The order of these scripts is intended to be introductory for learning, but it is not fixed. In real analyses, 
> you should move back and forth between steps, revisit earlier decisions, try different parameters.
> Don’t be afraid to experiment, explore the objects, and get your 
> hands dirty. That is often the most effective way to learn how these tools work.




