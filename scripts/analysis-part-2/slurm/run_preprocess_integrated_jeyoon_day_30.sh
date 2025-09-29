#!/bin/bash
#SBATCH --job-name=ventral_sosrs_preprocess_integrated_day30
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/preprocess_integrated_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/preprocess_integrated_%j.err
#SBATCH --partition=largemem
#SBATCH --time=6:00:00
#SBATCH --mem=1000G
#SBATCH --cpus-per-task=8

source ~/.bashrc
conda activate mindscape-env

export MINDSCAPE_INPUT_DIR="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_only_normalized_jeyoon_day_30"
export MINDSCAPE_OUTPUT_DIR="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_preprocessed_integrated_day_30"

Rscript /home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/mindscape_preprocess_integrated_jeyoon_day30.R
