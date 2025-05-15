#!/bin/bash

###############################################################################
# MindScape Interactive Environment Loader
#
# This script is meant to be used by anyone working interactively on the
# UMich Great Lakes HPC cluster, including in VS Code terminals or OnDemand.
#
# PURPOSE:
#   This script ensures that Conda is correctly initialized in your current
#   terminal and that the "mindscape-env" environment is activated.
#
# WHY THIS IS NEEDED:
#   When you load the Anaconda module, Conda becomes available, but it will
#   NOT fully work with the `conda activate` command unless we also source
#   the Conda shell initialization script.
#
#   Without this step, you'll get the error:
#       "CondaError: Run 'conda init' before 'conda activate'"
#
#   This script fixes that without needing to modify your .bashrc
#
# HOW TO USE:
#   Run this script every time you open a new terminal on Great Lakes:
#
#       source scripts/load_env_interactive.sh
#
#   NOTE: Do NOT use 'bash' to run this script, or the environment will not stay active.
###############################################################################

# Safety check: if this script is run with 'bash', it will not work as intended.
# We test if the script is being sourced (return works only when sourced).
(return 0 2>/dev/null)
if [ $? -ne 0 ]; then
  echo ""
  echo "WARNING: This script must be sourced, not executed with 'bash'."
  echo "Please run it like this:"
  echo ""
  echo "    source scripts/load_env_interactive.sh"
  echo ""
  echo "If you use 'bash', the Conda environment will NOT stay active after the script finishes."
  echo ""
  exit 1
fi

# Step 1: Load the UMich Anaconda module
module load python3.11-anaconda/2024.02

# Step 2: Initialize Conda for your shell session
source "$(conda info --base)/etc/profile.d/conda.sh"

# Step 3: Activate the mindscape Conda environment
conda activate mindscape-env

# Step 4: Confirm the environment is active
echo ""
echo "Environment 'mindscape-env' is now active and ready to use."
