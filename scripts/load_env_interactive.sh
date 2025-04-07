#!/bin/bash

###############################################################################
# MindScape Interactive Session Setup (Great Lakes + VS Code)
#
# This script is intended for interactive use only. It loads the necessary
# modules and activates the mindscape Conda environment in your current terminal.
#
# Usage (in terminal):
#   bash scripts/load_env_interactive.sh
###############################################################################

module load python3.11-anaconda/2024.02
conda activate mindscape-env
