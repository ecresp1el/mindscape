#!/bin/bash

# ==============================================================================
# SLURM Job Script for MindScape - ventral_sosrs_workflow_2025
# ==============================================================================
# HOW TO USE THIS SLURM SCRIPT (for UMich Great Lakes)
#
# This script runs 6 SLURM array jobs, each processing one sample using the
# MindScape R script for FLEX organoid scRNA-seq data. You must:
#
# 1. Place your R script at:
#    /nfs/turbo/umms-parent/Manny_test/MindScape/mindscape/workflows/ventral_sosrs_workflow_2025/mindscape_process_sample.R
#
# 2. Ensure Cell Ranger outputs exist at:
#    /nfs/turbo/umms-parent/Manny_test/10496-MW-reanalysis/outs/per_sample_outs/<sample_id>/count/
#
# 3. Submit the job using:
#    sbatch run_ventral_sosrs_slurm.sh
#
# Each task will:
# - Run the R script on one sample from SAMPLE_LIST
# - Save outputs to:
#     /nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/<sample_id>/
# - Save logs to:
#     logs/sosrs_<jobid>_<taskid>.out and .err
#
# You can monitor jobs with:
#     squeue -u $USER
#     tail -f logs/sosrs_<jobid>_<taskid>.out
# ==============================================================================
#
# This script runs per-sample processing for a list of single-cell RNA-seq
# datasets produced using Cell Ranger multi. Each sample is processed using an
# R script that expects the following environment variables:
#
#   - MINDSCAPE_SINGLE_SAMPLE: ID of the sample to process (e.g., "10496-MW-1")
#   - MINDSCAPE_OUTPUT_DIR: folder where results should be saved
#
# IMPORTANT: This script is meant to run as a SLURM array job where each task
# handles one sample independently. Logs and outputs are saved by sample.
#
# Update the following as needed:
#   - --account       → your UMich HPC account or grant
#   - SAMPLE_LIST     → list of sample IDs you wish to process
#   - SCRIPT_PATH     → path to your per-sample R processing script
# ==============================================================================

#SBATCH --job-name=ventral_sosrs                         # Job name shown in job queue
#SBATCH --account=parent0                               # UMich HPC billing account (confirmed)
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/sosrs_%A_%a.out                    # Log file (STDOUT)
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/sosrs_%A_%a.err                     # Log file (STDERR)
#SBATCH --time=8:00:00                                   # Time allocation per sample
#SBATCH --mem=32G                                        # Memory per task
#SBATCH --cpus-per-task=8                                # CPU cores per task
#SBATCH --array=0-5                                      # Number of samples (0-indexed)
#SBATCH --mail-type=FAIL                                 # Email only on failure
#SBATCH --mail-user=elcrespo@umich.edu                   # Change to your email

# ------------------------------------------------------------------------------
# Environment and software setup
# ------------------------------------------------------------------------------

# Initialize Conda (needed to use 'conda activate')
source ~/.bashrc

# Activate your fully configured MindScape Conda environment
conda activate mindscape-env

# ------------------------------------------------------------------------------
# Define base locations for code and output
# ------------------------------------------------------------------------------

BASE_DIR="/nfs/turbo/umms-parent/Manny_test"

# Folder containing your custom R script
SCRIPT_DIR="/home/elcrespo/Desktop/githubprojects/MindScape/mindscape/workflows/ventral_sosrs_workflow_2025"

# Path to the per-sample R script that handles processing logic
SCRIPT_PATH="$SCRIPT_DIR/mindscape_process_sample.R"

# Define a structured output location outside the Git repo, named after the job
OUTPUT_DIR="$BASE_DIR/ventral_sosrs_output"

# Log output folder now set to Turbo output location for consistency
LOG_DIR="$OUTPUT_DIR/logs"

mkdir -p $LOG_DIR $OUTPUT_DIR

# ------------------------------------------------------------------------------
# Define list of sample IDs — must match Cell Ranger per_sample_outs directory
# ------------------------------------------------------------------------------

SAMPLE_LIST=(10496-MW-1 10496-MW-2 10496-MW-3 10496-MW-4 10496-MW-5 10496-MW-6)

# Select sample for this SLURM array task
SAMPLE_ID=${SAMPLE_LIST[$SLURM_ARRAY_TASK_ID]}

# Export variables so R script knows what sample and output path to use
export MINDSCAPE_SINGLE_SAMPLE=$SAMPLE_ID
export MINDSCAPE_OUTPUT_DIR="$OUTPUT_DIR/$SAMPLE_ID"
mkdir -p $MINDSCAPE_OUTPUT_DIR

# ------------------------------------------------------------------------------
# Echo useful runtime info to log
# ------------------------------------------------------------------------------
echo "🧠 SLURM Job ID: $SLURM_JOB_ID"
echo "🧠 SLURM Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "🧠 Running sample: $SAMPLE_ID"
echo "🧠 Output directory: $MINDSCAPE_OUTPUT_DIR"
echo "🧠 Script being run: $SCRIPT_PATH"
echo "🧠 Conda environment: $CONDA_DEFAULT_ENV"

# ------------------------------------------------------------------------------
# Run the R script for the selected sample
# ------------------------------------------------------------------------------

Rscript $SCRIPT_PATH