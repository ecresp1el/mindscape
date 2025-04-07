#!/bin/bash

###############################################################################
# MindScape Environment Setup Script (Local or Non-UMich)
#
# This version assumes you are running on a machine where you have conda
# installed directly (not via module system).
#
# Usage:
#   bash scripts/setup_mindscape_env.sh
###############################################################################

# Step 1: Define environment name and YAML file path
ENV_NAME="mindscape-env"
ENV_FILE="conda-environments/MINDSCAPE.yml"

# Step 2: Check that conda is available
if ! command -v conda &> /dev/null; then
  echo "ERROR: Conda is not installed or not on your PATH."
  echo "Please install Anaconda or Miniconda before running this script."
  exit 1
fi

# Step 3: Check for the YAML file
if [ ! -f "$ENV_FILE" ]; then
  echo "ERROR: Cannot find environment file: $ENV_FILE"
  exit 1
fi

# Step 4: Remove existing environment (if any)
if conda info --envs | grep -q "$ENV_NAME"; then
  echo "Environment '$ENV_NAME' already exists. Removing it..."
  conda remove --name $ENV_NAME --all -y
fi

# Step 5: Create the environment
echo "Creating conda environment from $ENV_FILE..."
conda env create -f $ENV_FILE

# Step 6: Activate and test the environment
echo "Activating and testing environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME
python -c "import scanpy; print('Python environment setup verified.')"

echo "MindScape environment '$ENV_NAME' is ready."
