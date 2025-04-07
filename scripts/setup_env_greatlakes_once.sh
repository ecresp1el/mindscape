#!/bin/bash

###############################################################################
# MindScape Launch Script
#
# This script ensures all required modules and environments are loaded before
# running any downstream analysis. It is fully reproducible and compatible
# with UMich's Great Lakes HPC.
#
# Usage:
#   bash scripts/launch_env_with_modules.sh
###############################################################################

# Step 1: Load required UMich modules
echo "Loading UMich Great Lakes modules..."
module load python3.11-anaconda/2024.02
module load Bioinformatics
module load cellranger/9.0.1

# Step 2: Activate the MindScape conda environment
ENV_NAME="mindscape-env"
echo "Activating conda environment: $ENV_NAME"
source activate $ENV_NAME

# Step 3: Verify Python environment works
echo "Verifying Python environment..."
python -c "import scanpy; print('Python packages loaded successfully.')"

# Step 4: Verify Cell Ranger CLI is available
echo "Verifying Cell Ranger..."
cellranger --version || {
  echo "ERROR: Cell Ranger is not available in this shell."
  exit 1
}

# Step 5: Placeholder for your analysis command
# You can run anything you want after this point, like:
# python scripts/run_analysis.py
# or:
# cellranger count --id=sample1 ...

echo "Environment and modules are ready. Add your workflow commands here."
