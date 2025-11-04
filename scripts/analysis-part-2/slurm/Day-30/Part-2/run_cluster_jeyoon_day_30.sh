#!/bin/bash
#SBATCH --job-name=ventral_sosrs_cluster_day30
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/cluster_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/cluster_%j.err
#SBATCH --partition=standard
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8

source ~/.bashrc
conda activate mindscape-env

export MINDSCAPE_PREPROCESSED_RDS="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_preprocessed_integrated_day_30/preprocessed_integrated_jeyoon_day_30.rds"
export MINDSCAPE_OUTPUT_DIR="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_clustered_day_30"

Rscript /home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/mindscape_cluster_jeyoon_day30.R
