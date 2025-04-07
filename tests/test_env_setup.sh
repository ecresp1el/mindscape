#!/bin/bash

###############################################################################
# MindScape Environment Test Script
#
# This script checks your system to help you determine whether you're ready
# to run MindScape's setup scripts. It works for both local machines and
# the UMich Great Lakes HPC cluster.
#
# Usage:
#   bash scripts/test_env_setup.sh
###############################################################################

echo "==============================================="
echo "MindScape Environment Setup Test - Starting"
echo "==============================================="
echo ""

# Step 1: Detect host system
if hostname | grep -q "gl"; then
  echo "Detected system: UMich Great Lakes HPC"
  ON_GREAT_LAKES=true
else
  echo "Detected system: Local or non-UMich machine"
  ON_GREAT_LAKES=false
fi

echo ""

# Step 2: Check if conda is available
if command -v conda &> /dev/null; then
  echo "Conda is available in this terminal session."
  CONDA_OK=true
else
  echo "Conda is NOT currently available."
  CONDA_OK=false
fi

# Step 3: If on Great Lakes, check if Anaconda module is available
if $ON_GREAT_LAKES; then
  echo ""
  echo "Checking for Anaconda module on Great Lakes..."
  if module avail python3.11-anaconda/2024.02 2>&1 | grep -q "python3.11-anaconda"; then
    echo "Anaconda module is available and can be loaded."
    MODULE_OK=true
  else
    echo "Anaconda module is NOT available. You may not be on a compute or login node."
    MODULE_OK=false
  fi
fi

# Step 4: Check if the mindscape environment already exists
echo ""
if command -v conda &> /dev/null && conda info --envs | grep -q "mindscape-env"; then
  echo "Conda environment 'mindscape-env' already exists."
  ENV_EXISTS=true
else
  echo "Conda environment 'mindscape-env' does NOT exist."
  ENV_EXISTS=false
fi

# Step 5: Final recommendation
echo ""
echo "Summary and Suggested Next Step:"
if ! $CONDA_OK && $ON_GREAT_LAKES && $MODULE_OK && ! $ENV_EXISTS; then
  echo "- Run the setup script for Great Lakes:"
  echo "    bash scripts/setup_env_greatlakes_once.sh"
elif ! $CONDA_OK && $ON_GREAT_LAKES && $MODULE_OK && $ENV_EXISTS; then
  echo "- Conda is not currently available, but your environment exists."
  echo "  Run the interactive script to enter your environment:"
  echo "    bash scripts/load_env_interactive.sh"
elif ! $CONDA_OK && ! $ON_GREAT_LAKES; then
  echo "- Install Anaconda or Miniconda on your local machine and rerun this test."
elif $CONDA_OK && ! $ENV_EXISTS; then
  echo "- Conda is available. Run the local setup script:"
  echo "    bash scripts/setup_env_local.sh"
elif $CONDA_OK && $ENV_EXISTS; then
  echo "- Conda and the environment are both ready. Run:"
  echo "    conda activate mindscape-env"
else
  echo "- Please check the messages above and troubleshoot as needed."
fi

