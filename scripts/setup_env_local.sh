#!/bin/bash

###############################################################################
# MindScape Environment Setup Script (Local or Non-UMich)
#
# Creates a clean Conda environment from the project's YAML file.
#
# Designed for use on:
#   - Personal machines (with Miniconda or Anaconda installed)
#   - HPC clusters like UMich Great Lakes (tries to auto-load Anaconda module)
#
# Usage:
#   bash scripts/setup_env_local.sh
###############################################################################

ENV_NAME="mindscape-env"
ENV_FILE="conda-environments/MINDSCAPE.yml"

# Step 1: Ensure conda is available
if ! command -v conda &> /dev/null; then
  echo "WARNING: Conda command not found."
  echo "Trying to load the Anaconda module (this may work on UMich Great Lakes)..."

  module load python3.11-anaconda/2024.02

  if ! command -v conda &> /dev/null; then
    echo ""
    echo "ERROR: Conda is still not available."
    echo "If you're on Great Lakes, run: module load python3.11-anaconda/2024.02"
    echo "If you're on your own computer, install Anaconda or Miniconda."
    exit 1
  fi

  echo "Conda successfully loaded from the module system."
else
  echo "Conda is available."
fi

# Step 2: Ensure the YAML file exists
if [ ! -f "$ENV_FILE" ]; then
  echo "ERROR: Environment file not found: $ENV_FILE"
  echo "Are you running this from the MindScape project root?"
  exit 1
fi

# Step 3: Remove existing environment (if any)
if conda info --envs | grep -q "$ENV_NAME"; then
  echo "Removing existing environment: $ENV_NAME"
  conda remove --name $ENV_NAME --all -y
fi

# Step 4: Create the environment
echo "Creating new Conda environment: $ENV_NAME"
conda env create -f $ENV_FILE

# Step 5: Activate and test environment
echo "Activating environment and verifying installation..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME

if python -c "import scanpy" &> /dev/null; then
  echo "Environment setup complete. scanpy is available."
else
  echo "Environment created, but scanpy could not be imported."
fi
