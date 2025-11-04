#!/bin/bash
# ==============================================================================
# MindScape Environment Setup Script for Great Lakes HPC (UMich)
# ==============================================================================

# Exit the script immediately if any command fails
set -e

echo "Starting MindScape environment setup..."

# ------------------------------------------------------------------------------
# Step 1: Define install location and detect Git root
# ------------------------------------------------------------------------------
# Set the directory where Miniconda will be installed (in your home directory)
INSTALL_DIR="$HOME/miniconda3"
# The name of the environment YAML file (relative to repo root)
ENV_YML_NAME="conda-environments/MINDSCAPE.yml"  # Relative path from repo root

# Find the top-level directory of the current Git repository
REPO_ROOT=$(git rev-parse --show-toplevel)
# Full path to the environment YAML file
YML_PATH="${REPO_ROOT}/${ENV_YML_NAME}"

# ------------------------------------------------------------------------------
# Step 2: Install Miniconda if not already installed
# ------------------------------------------------------------------------------
# Check if Miniconda is already installed in the desired location
if [ ! -d "$INSTALL_DIR" ]; then
    echo " Installing Miniconda at $INSTALL_DIR"
    # Download the Miniconda installer script
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    # Run the installer in batch mode (no prompts), install to INSTALL_DIR
    bash miniconda.sh -b -p "$INSTALL_DIR"
    # Remove the installer script after installation
    rm miniconda.sh
else
    echo "Miniconda already installed at $INSTALL_DIR"
fi

# ------------------------------------------------------------------------------
# Step 3: Initialize Conda and update shell config
# ------------------------------------------------------------------------------
echo "Initializing Conda for bash"
# Source the conda.sh script to enable 'conda' command in this shell
source "$INSTALL_DIR/etc/profile.d/conda.sh"

# ------------------------------------------------------------------------------
# Step 3b: Ensure mamba is available
# ------------------------------------------------------------------------------
if ! command -v mamba &> /dev/null; then
    echo "🔧 Mamba not found. Creating temporary bootstrap environment with mamba..."
    conda create -n mamba-bootstrap -y -c conda-forge mamba
    conda activate mamba-bootstrap
else
    echo "Mamba is already available."
fi

# Only initialize Conda in bash if it hasn't been done before
if ! grep -q "conda initialize" ~/.bashrc; then
    # Adds Conda initialization code to your ~/.bashrc
    conda init bash
fi

# Disable automatic activation of the base environment (recommended by UMich HPC)
conda config --set auto_activate_base false

# ------------------------------------------------------------------------------
# Step 4: Configure Conda channels for bioinformatics work
# ------------------------------------------------------------------------------
echo "Configuring conda channels"
# Add conda-forge channel for up-to-date packages
conda config --add channels conda-forge
# Add bioconda channel for bioinformatics software
conda config --add channels bioconda
# Enforce strict channel priority for reproducibility
conda config --set channel_priority strict

# ------------------------------------------------------------------------------
# Step 5: Create the environment from MINDSCAPE.yml
# ------------------------------------------------------------------------------
# Check if the environment YAML file exists
if [ -f "$YML_PATH" ]; then
    echo "Found environment file at $YML_PATH"
    echo "Creating environment from YAML..."

    # Check if the environment 'mindscape-env' exists, update or create accordingly
    if [ -d "$HOME/miniconda3/envs/mindscape-env" ]; then
        echo "Environment 'mindscape-env' already exists. Updating it from YAML..."
        mamba env update --name mindscape-env --prune -f "$YML_PATH"
    else
        echo "Creating new environment 'mindscape-env' from YAML..."
        mamba env create -f "$YML_PATH"
    fi
else
    echo "ERROR: Could not find $YML_PATH"
    exit 1
fi


# ------------------------------------------------------------------------------
# Step 6: Install R packages from GitHub (e.g., SeuratDisk)
# ------------------------------------------------------------------------------
echo "Installing R GitHub packages (e.g., SeuratDisk)..."
conda activate mindscape-env
Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes"); remotes::install_github("mojaveazure/seurat-disk")'

echo "✅ Setup complete. Run 'conda activate mindscape-env' to begin."