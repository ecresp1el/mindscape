#!/usr/bin/env bash
# run_pipeline.sh
# Description: Robust wrapper to launch Snakemake from this GitHub repo using a shared Turbo project directory.

set -euo pipefail

echo "🚀 Launching MindScape pipeline setup..."

#############################################
# STEP 1: Run first rule: create_project
#############################################
echo "📦 Creating project structure..."

# Use config file from this GitHub repo
GITHUB_CONFIG="config/config.yaml"

snakemake --snakefile workflow/Snakefile \
          --configfile "$GITHUB_CONFIG" \
          --cores 1 \
          --printshellcmds \
          --rerun-incomplete \
          create_project

#############################################
# STEP 2: Read updated project_path from config
#############################################
PROJECT_PATH=$(python -c "import yaml; print(yaml.safe_load(open('$GITHUB_CONFIG'))['project_path'])")
UPDATED_CONFIG="$PROJECT_PATH/config/config.yaml"

echo "✅ Project directory created at: $PROJECT_PATH"

#############################################
# STEP 3: Run full pipeline using updated config
#############################################
echo "🔁 Running remaining workflow steps..."

snakemake --snakefile /home/elcrespo/Desktop/githubprojects/MindScape/mindscape_snakemake/workflow/Snakefile \
          --configfile "$UPDATED_CONFIG" \
          --directory "$PROJECT_PATH" \
          --cores 1 \
          --printshellcmds \
          --rerun-incomplete