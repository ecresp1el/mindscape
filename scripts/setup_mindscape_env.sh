#!/bin/bash

###############################################################################
# MindScape Environment Setup Script
#
# This script sets up the 'mindscape-env' Conda environment using a YAML file
# and loads the appropriate UMich HPC modules, including Cell Ranger.
#
# Requirements:
# - You are logged into the Great Lakes cluster.
# - You have cloned the MindScape GitHub repository.
#
# Usage:
#   bash scripts/setup_mindscape_env.sh
###############################################################################

# Step 1: Load modules provided by UMich ITS on the Great Lakes cluster
# This ensures we are using UM-maintained software installations
module load python3.11-anaconda/2024.02
module load cellranger/9.0.1

# Step 2: Define environment details
ENV_NAME="mindscape-env"
ENV_FILE="conda-environments/MINDSCAPE.yml"

# Step 3: Ensure the YAML file exists
if [ ! -f "$ENV_FILE" ]; then
  echo "ERROR: Cannot find environment file: $ENV_FILE"
  echo "Make sure you are in the MindScape project root and cloned from GitHub."
  exit 1
fi

# Step 4: Remove existing environment (if any) to ensure a clean setup
if conda info --envs | grep -q "$ENV_NAME"; then
  echo "Environment '$ENV_NAME' already exists. Removing it..."
  conda remove --name $ENV_NAME --all -y
fi

# Step 5: Create the Conda environment from the provided YAML file
echo "Creating conda environment '$ENV_NAME' from $ENV_FILE..."
conda env create -f $ENV_FILE

# Step 6: Activate and test the environment
echo "Activating and testing environment..."
source activate $ENV_NAME
python -c "import scanpy; print('Python environment setup verified.')"

# Step 7: Verify Cell Ranger is loaded and usable
echo "Verifying Cell Ranger installation..."
cellranger --version || {
  echo "ERROR: Cell Ranger is not available. Check that the module loaded correctly."
  exit 1
}

echo "Environment setup complete. Both Python and Cell Ranger are ready to use."
