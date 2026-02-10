#############SECTION 4: HARMONY INTEGRATION####################################################

# --- Merge vs. Integration: When to Use Which ---
# Merge (Seurat::merge):
#   • Combines multiple Seurat objects into one.
#   • Ideal when samples share similar technical conditions and batch effects are minor.
# Integration (Harmony):
#   • Actively corrects technical or batch differences between datasets.
#   • Essential when combining data from different experimental dates, sample preps, timepoints.
#   • Aligns shared biological signals across batches, improving clustering and visualization.
#   • For some datasets, it can wash away TRUE biological differences (e.g. integration healthy vs disease conditions).
# Decide what is best for your dataset!!


#Start by using your original merged object from section 3 
#ensure the merged obj did not go through any manipulation
int_obj <- readRDS("merged_raw_obj.rds") #"merged_raw_obj.rds"


#Normalize using SCTransform 
int_obj<- SCTransform(merged_obj,assay = "Xenium", verbose = FALSE)

#Run PCA 
int_obj<- RunPCA(int_obj,features = rownames(merged_obj), npcs = 100,
                 verbose = FALSE)

#Select PCs via ElbowPlot (or other methods discussed above)
ElbowPlot(int_obj,ndims = 100) 

# Integrate samples with Harmony
#    - orig.reduction: PCA embeddings to correct
#    - new.reduction: stores the batch-corrected embeddings
#    - normalization.method: use SCT-normalized data
int_obj<- IntegrateLayers(
  object = int_obj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  normalization.method = "SCT",
  verbose = TRUE)

#Run UMAP, find neighbours and clustering
int_obj<- RunUMAP(int_obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
int_obj<- FindNeighbors(int_obj, reduction = "harmony", dims = 1:20)
int_obj<- FindClusters(int_obj, resolution = 0.1)