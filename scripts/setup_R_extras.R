

# ------------------------------------------------------------------------------
# setup_R_extras.R
#
# This script installs additional R packages from GitHub that are required
# for cross-language functionality in MindScape, specifically SeuratDisk.
# SeuratDisk enables conversion between Seurat (R) and Scanpy (Python)
# compatible formats such as .h5Seurat and .h5ad, which is essential for
# workflows that span both environments.
# ------------------------------------------------------------------------------

# Ensure the 'remotes' package is available. It allows us to install R packages from GitHub.
# If it is not installed, install it from CRAN.
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Use remotes to install SeuratDisk from its official GitHub repository.
# This is necessary because SeuratDisk is not available on CRAN or Conda,
# but it is required for .h5Seurat ↔ .h5ad conversion.
remotes::install_github("mojaveazure/seurat-disk")

# Confirm the installation
cat("✅ SeuratDisk successfully installed.\n")