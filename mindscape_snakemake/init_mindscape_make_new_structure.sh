#!/bin/bash

# Exit if any command fails
set -e

ROOT_DIR=$(pwd)

echo "📁 Initializing Snakemake project structure in: $ROOT_DIR"

mkdir -p workflow/rules
mkdir -p workflow/envs
mkdir -p workflow/scripts
mkdir -p workflow/notebooks
mkdir -p workflow/report
mkdir -p config
mkdir -p results
mkdir -p resources

touch README.md
touch LICENSE.md
touch workflow/Snakefile
touch workflow/rules/module1.smk
touch workflow/rules/module2.smk
touch workflow/envs/tool1.yaml
touch workflow/envs/tool2.yaml
touch workflow/scripts/script1.py
touch workflow/scripts/script2.R
touch workflow/notebooks/notebook1.py.ipynb
touch workflow/notebooks/notebook2.r.ipynb
touch workflow/report/plot1.rst
touch workflow/report/plot2.rst
touch config/config.yaml
touch config/some-sheet.tsv

echo "✅ Snakemake project scaffold created at $ROOT_DIR"