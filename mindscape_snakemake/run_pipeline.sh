# run_pipeline.sh
# Description: Robust wrapper to launch Snakemake from this GitHub repo using a shared Turbo project directory

set -euo pipefail

# Step 1: Inform user of start
echo "🚀 Launching MindScape pipeline setup..."

# Step 2: Run the first Snakemake rule to create the project directory
echo "📦 Creating project structure..."
snakemake --snakefile workflow/Snakefile \
          --configfile config/config.yaml \
          --cores 1 \
          --printshellcmds \
          --rerun-incomplete \
          create_project

# Step 3: Extract project path from the config file for user reference
CONFIG_FILE="config/config.yaml"
PROJECT_PATH=$(python -c "import yaml; print(yaml.safe_load(open('$CONFIG_FILE'))['project_path'])")

# Step 4: Inform the user
echo "✅ Project directory created at: $PROJECT_PATH"

# Step 5: Run the rest of the workflow from this GitHub location
echo "🔁 Running remaining workflow steps..."
snakemake --snakefile workflow/Snakefile \
          --configfile config/config.yaml \
          --cores 1 \
          --printshellcmds \
          --rerun-incomplete
