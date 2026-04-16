# =========================================================================
# ## Load Libraries
#
# First 
# - install.packages('BiocManager')
# 
# as it is a wrapper for both:
#
# - install.packages()
# - remotes::install_github()
# 
# So, it works as a "best of all worlds" approach for bioinformatics packages.
#
#
# Before running, in Rstudio Menu bar follow:
#
# - Session -> Restart R 
# - Or use shortcut: (Ctrl + Shift + F10)
#
# To avoid conflicts due packages already loaded.
#
# You can now run the whole script using shortcut (Ctrl + Shift + Enter)
# =========================================================================


# Define the complete list of required packages
packages <- c(
  # Data Import & Organization
  "here", "data.table", "arrow",
  
  # Data Manipulation
  "dplyr", "tidyverse",
  
  # Spatial & Single-Cell Processing
  "Seurat", "spacexr", "glmGamPoi",
  
  # Integration & Batch Correction
  "harmony",
  
  # Visualization
  "ggplot2", "patchwork", "RColorBrewer", "scCustomize", 
  "scattermore", "paletteer", "Polychrome",
  
  # Parallelization & Performance
  "future",
  
  # SingleR Analysis
  "SingleR", "rlang", "geometry", "purrr", "DESeq2"
)

# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Function to check, install, and update packages
setup_packages_only <- function(pkgs) {
  message("Starting package installation and update process...")
  
  # Install if package is missing
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing missing package:", pkg))
      # BiocManager::install for both CRAN and Bioconductor packages
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      message(paste("Package already installed:", pkg))
    }
  }
  
  # Update all specified packages to their latest versions
  message("Checking for package updates...")
  BiocManager::install(pkgs, update = TRUE, ask = FALSE)
  
  message("================================\n", 
          "Installation and update process complete!\n", 
          "================================")
}

# Run the installation and update function
setup_packages_only(packages)

# remove function and list from memory for efficiency
rm(packages, setup_packages_only)

gc()
