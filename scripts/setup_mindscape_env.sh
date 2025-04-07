#!/bin/bash

###############################################################################
# MindScape Environment Setup Script
#
# Sets up the mindscape-env Conda environment on UMich Great Lakes HPC
# after cloning the MindScape repository from GitHub.
#
# Assumptions:
# - You are running this script from within a cloned GitHub repository.
# - The conda-environments/MINDSCAPE.yml file exists.
#
# How to use:
#   git clone https://github.com/ecresp1el/MindScape.git
#   cd MindScape
#   bash scripts/setup_mindscape_env.sh
###############################################################################

# Step 1: Load UMich-provided Anaconda Python compiler module
module load python3.11-anaconda/2024.02

# Step 2: Define environment name and YAML file path
ENV_NAME="mindscape-env"
ENV_FILE="conda-environments/MINDSCAPE.yml"

# Step 3: Check that the YAML file exists
if [ ! -f "$ENV_FILE" ]; then
  echo "Error: Cannot find environment file: $ENV_FILE"
  echo "Make sure you are in the root of the cloned GitHub repository."
  exit 1
fi

# Step 4: If the environment already exists, remove it
if conda info --envs | grep -q "$ENV_NAME"; then
  echo "Environment '$ENV_NAME' already exists. Removing it..."
  conda remove --name $ENV_NAME --all -y
fi

# Step 5: Create the environment from the YAML file
echo "Creating conda environment from $ENV_FILE..."
conda env create -f $ENV_FILE

# Step 6: Test that key packages are installed correctly
echo "Activating and testing environment..."
source activate $ENV_NAME
python -c "import scanpy, cellwrangler; print('Environment setup verified successfully.')"

echo "Environment '$ENV_NAME' is ready."
