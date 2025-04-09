#!/bin/bash

###############################################################################
# Download and Set Up Human GRCh38 Reference for Cell Ranger
#
# This script ensures that the required 10x Genomics reference is available
# in a reproducible location on Turbo storage.
#
# It will:
# - Create /nfs/turbo/umms-parent/Manny_human_ref if it doesn't exist
# - Download and extract the GRCh38 2020-A prebuilt reference from 10x Genomics
# - Print progress messages
#
# HOW TO RUN:
#   bash scripts/setup_human_reference.sh
###############################################################################

REF_PARENT="/nfs/turbo/umms-parent"
REF_DIR="${REF_PARENT}/Manny_human_ref"
REF_SUBDIR="${REF_DIR}/refdata-gex-GRCh38-2020-A"
TARBALL="refdata-gex-GRCh38-2020-A.tar.gz"
DOWNLOAD_URL="https://cf.10xgenomics.com/supp/cell-exp/${TARBALL}"

echo "Checking if reference directory already exists..."

if [ -d "$REF_SUBDIR" ]; then
  echo "Reference already exists at: $REF_SUBDIR"
else
  echo "Creating reference folder: $REF_DIR"
  mkdir -p "$REF_DIR"
  cd "$REF_DIR" || { echo "Error: Failed to cd into $REF_DIR"; exit 1; }

  echo "Downloading reference from 10x Genomics..."
  wget "$DOWNLOAD_URL"

  echo "Extracting reference files..."
  tar -xzvf "$TARBALL"

  echo "Cleaning up downloaded archive..."
  rm "$TARBALL"

  echo "Reference successfully set up at: $REF_SUBDIR"
fi
