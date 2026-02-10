# PCA


Principal Component Analysis (PCA) is a dimensionality reduction technique 
used in data science to simplify large datasets by transforming numerous, 
potentially correlated variables into a smaller set of uncorrelated variables 
called "principal components," while retaining most of the original information.
It is highly effective for data visualization, noise reduction, and speeding up 
 algorithms. 
 

Most Relevant arguments:


- npcs: No. of Principal Components you want
- features: Features to compute PCA on. By default uses variable features for the Assay 
(We use VariableFeatures(MAID_obj) for SCT variable features)
 
Link to Seurat's RunPCA implementation: [here](https://satijalab.org/seurat/reference/runpca)