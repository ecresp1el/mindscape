#!/bin/bash

# run_pipeline.sh
# Description: Wrapper script to launch Snakemake in the newly created project directory on Turbo

# Exit on any error
set -euo pipefail

# Step 1: Run Snakemake to create the project directory
snakemake -j 1 -p create_project

# Step 2: Parse the project_path from the config file
CONFIG_FILE="config/config.yaml"
PROJECT_PATH=$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['project_path'])")

# Step 3: Echo the detected project path
echo "📁 Switching to project directory: $PROJECT_PATH"

# Step 4: Change to that directory
cd "$PROJECT_PATH"

# Step 5: Re-run Snakemake from inside the new directory
snakemake -j 1 -p
