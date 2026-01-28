---
title: "Xenium Demo"
date: "01/28/26"
---

```
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width = 8,
  fig.height = 6
)
```

# Overview

TODO: The goal is not just to run code, but to **explain why each step exists** and what decisions users should think about during analysis.

---

# Prerequisites

- R (>= 4.2 recommended), Rstudio
- Load objects provided by UBC (`.rds`)
- 8GB of memory,(recommended memory >=16GB)

```
{r libraries}
library(here)
library(data.table)
library(arrow)

library(dplyr)
library(tidyverse)

library(Seurat)
library(spacexr)
library(glmGamPoi)
library(harmony)

library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scCustomize)
library(scattermore)
library(Polychrome)

library(future)
library(SingleR)
library(rlang)
library(geometry)
library(purrr)

options(future.globals.maxSize = 250000 * 1024^2)
```

If a library does not load, try this:

- install.packages("_library_name_")

If that gives an error(eg for spacexr, and glmGamPoi), install BioCManager by running:

  - install.packages('BiocManager')
  - BiocManager::install('_library_name_')

---

## Load Xenium Seurat Object

> **Important naming note**: Avoid naming your R objects the same as Xenium FOV identifiers or `orig.ident` values.

```
setwd("C:/_Working_Directory_") # Eg. D:/work/Xenium/Input_RDS 

MAID_obj <- readRDS("M_70_Seurat_obj.rds")
MAID_obj <- UpdateSeuratObject(MAID_obj)
```

---

# PRELIMINARY ANALYSIS:

---


## Quality Control (QC)

We start by examining transcript depth and feature complexity per cell.

```
# nFeature_Xenium
ggplot(MAID_obj@meta.data, aes(x = nFeature_Xenium)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  scale_y_log10() +
  theme_classic() +
  labs(title = "nFeature_Xenium",
       x = "Number of detected features",
       y = "Cell count")

# nCount_Xenium
ggplot(MAID_obj@meta.data, aes(x = nCount_Xenium)) +
  geom_histogram(bins = 25, fill = "firebrick", color = "black") +
  scale_y_log10() +
  theme_classic() +
  labs(title = "nCount_Xenium",
       x = "Total transcript counts",
       y = "Cell count")
```

Filtering is dataset-dependent and should be guided by these distributions.

```
# Filter cells based on QC thresholds
MAID_obj <- subset(MAID_obj, subset = nFeature_Xenium > 2 & nCount_Xenium > 0) 
# Note: to compare multiple QC filters for downstream analysis
# assign filtered results to new object names 
# (e.g., MAID_obj_filtered) instead of overwriting MAID_obj.
# MAID_obj_filtered <- subset(MAID_obj, subset = nFeature_Xenium > 0 & nCount_Xenium > 0)
```

---

## Normalization with SCTransform

SCTransform performs:

- Normalization
- Variance stabilization
- Feature selection

Link to documentation on SCTransform: https://satijalab.org/seurat/articles/sctransform_vignette.html
Link to paper on SCTransform: https://link.springer.com/article/10.1186/s13059-021-02584-9

```
MAID_obj <- SCTransform(MAID_obj, assay = "Xenium")
```

You can analyze the normalisation, using the following code, we're looking for more/less a flat line.
```

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
gene_attr_maid <- SCTResults(MAID_obj, slot = "feature.attributes", assay = "SCT")

residualVarPlot(gene_attr_maid, max_resvar = 10, pt_size = 4)

```

---

## PCA and Dimensionality Reduction

PCA reduces high-dimensional data to principal components

- features: use VariableFeatures from SCTransform
- npcs: number of PCs to compute (e.g., 50)

Balance is Key: 

- Including more PCs captures more variation in the data, but may introduce noise.
- Including fewer PCs reduces noise but may risk missing important biological signals.

```
MAID_obj <- RunPCA(
  MAID_obj,
  features = VariableFeatures(MAID_obj),
  npcs = 50,
  verbose = FALSE
)
```

Examine PCA using these:
```
ElbowPlot(MAID_obj, ndims = 50)
```

> **Discussion point for the demo:**

- Does the elbow flatten out?
- What does this imply about variance captured by later PCs?
- What are the consequences (biological and computational) of choosing too many or too few PCs?

```
DimHeatmap(MAID_obj, dims = 1:12, cells = 500, balanced = TRUE) 
```
Question the implementation
---

## UMAP and Clustering

Umap: Uniform Manifold Approximation and Projection (UMAP) is a modern 
non-linear dimensionality reduction technique used to project high-dimensional 
data into a lower-dimensional space.

```
- dims: PCs to include (e.g., 1:20)
MAID_obj <- RunUMAP(MAID_obj, dims = 1:24)

MAID_obj <- FindNeighbors(MAID_obj, dims = 1:24)

# - resolution: higher => more clusters
MAID_obj <- FindClusters(MAID_obj, resolution = 0.5)

DimPlot(
  MAID_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  cols = "polychrome"
) + DarkTheme() + coord_fixed()
```

---

## Visualization

There are many ways of visualizing expression data. A few examples are shown:

```
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
```

---

## Spatial Visualization

```
ImageFeaturePlot(MAID_obj, fov = "M_70", features = "SLC17A7")
```

---

## Subsetting Strategies

There are several ways to subset (dependent on dataset), including:

1. Cluster identity
2. Marker expression thresholds
3. Cropping

```{r subset}
#1
obj_sub <- subset(MAID_obj, idents = c(0, 2, 4, 6))
#2
obj_sub <- subset(MAID_obj, subset = SLC17A7 > 1)
#3
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
```

You can Visualize the new data using:


```
ImageDimPlot(
  obj_sub,
  fov = "M_70",
  cols = "polychrome",
  size = 0.5,
  dark.background = TRUE
) + DarkTheme()
```

---


# Summary

TODO



