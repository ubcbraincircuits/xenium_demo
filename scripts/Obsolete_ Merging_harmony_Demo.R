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

####### LOAD IN RDS OBJECTS

# Set the working directory to the folder containing your saved RDS objects
setwd("D:/work/Github_demo/xenium_demo/Data")

