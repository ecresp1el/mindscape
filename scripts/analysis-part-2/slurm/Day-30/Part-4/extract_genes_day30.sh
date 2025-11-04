#!/bin/bash

# ==============================================================================
# SLURM Job Script for MindScape - Extract Multiple Genes with Cell Identity
# ==============================================================================

#SBATCH --job-name=extract_genes_day30
#SBATCH --account=parent0
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/extract_genes_day30_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/extract_genes_day30_%j.err
#SBATCH --partition=standard
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu

# ------------------------------------------------------------------------------
# Load environment
# ------------------------------------------------------------------------------
source ~/.bashrc
conda activate mindscape-env

# ------------------------------------------------------------------------------
# Set paths (Day 30)
# ------------------------------------------------------------------------------
BASE_DIR="/nfs/turbo/umms-parent/Manny_test"

# Input CSV (DEGs_all_clusters.csv)
export DEG_INPUT_CSV="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_day_30/DEGs_all_clusters_res_0.2/DEGs_all_clusters.csv"

# Output directory
export DEG_OUTPUT="$BASE_DIR/ventral_sosrs_output_integrated_analysis_cell_cycle_reg_jeyoon_day_30/gene_extracts_res_0.2"
mkdir -p "$DEG_OUTPUT"

# Target genes + identities (order matters: gene[i] → identity[i])
export TARGET_GENES="SOX2,VIM,MKI67,ESCO2,DLX5,DLL1,DLX2,ASCL1,HES6,DLX5,DCX,TUBB3,\
NEUROD1,ERBB4,CXCR4,DCX,ARX,ADAMTS5,BEND4,DLGAP,KIT1,RAI2,CALB1,CUX2,ERBB4,LMO3,\
MAF,MAFB,NPY,SOX6,SST,ZEB2,PVALB,CALB2,GAD2,DLX5,DLX2,ISL1,DLX1,SLC32A1,SP9,DLX6,\
DLX6-AS1,KCNMB2,ESRRG,GAD1,MIR7-3HG,SIX3,POU3F4,DCBLD2,POU3F2,RGS2,NANOS1,ARX,ZFHX3,VGF,\
XPR1,SCG3,PAX6,GNG3,CD24,LHX1,GAP43,ECEL1,PCSK1N,STMN2,RBP1,PRKACB,SST,MEIS2,NEUROD6,TLE4,\
TBR2,FOXP2,TIS21,TBR1,FEZF2,SEMA3E,ADCYAP1,LHX2,CUX1,PLXND1,VIM,NES,HES1,MOXD1,TCIM,HES5,\
FAM107A,HOPX,TNC,EGFR,DLL3,PRRX1,AGMO,ZNF703,PDGFRA,ENPP6,MAG,RITZ,RIT2,BMPER,OLIG1,\
AQP4,GJA1,SB100B,GFAP,NKX2-1,TOX3,PLS3,MEF2C,MAF,ARX,ZEB2,GBX,ISL1,NKX2-1,LHX8,DLX,\
ZEB2,PLS3,ERBB4,CXCR4,LHX8,ISL1,ZIC1,GBX2,ETV1,NFIA,NKX2-1,ST18,ETV1,ANGPT2,MAFB,MEF2C"

export TARGET_IDENTITIES="Dividing Cells/Progenitors,Dividing Cells/Progenitors,Dividing Cells/Progenitors,Dividing Cells/Progenitors,Dividing Cells/Progenitors,\
Intermediate Progenitors,Intermediate Progenitors,Intermediate Progenitors,Intermediate Progenitors,Intermediate Progenitors,\
Newborn Neurons,Newborn Neurons,Newborn Neurons,\
Migratory Interneurons,Migratory Interneurons,Migratory Interneurons,Migratory Interneurons,\
Immature Cortical Interneurons,Immature Cortical Interneurons,Immature Cortical Interneurons,Immature Cortical Interneurons,Immature Cortical Interneurons,\
Cortical Interneurons,Cortical Interneurons,Cortical Interneurons,Cortical Interneurons,Cortical Interneurons,Cortical Interneurons,Cortical Interneurons,Cortical Interneurons,Cortical Interneurons,Cortical Interneurons,Cortical Interneurons,Cortical Interneurons,\
GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,\
GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,\
GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,\
GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,GABAergic Neurons,\
Excitatory Neurons/Progenitors,Excitatory Neurons/Progenitors,Excitatory Neurons/Progenitors,Excitatory Neurons/Progenitors,Excitatory Neurons/Progenitors,Excitatory Neurons/Progenitors,Excitatory Neurons/Progenitors,Excitatory Neurons/Progenitors,Excitatory Neurons/Progenitors,Excitatory Neurons/Progenitors,Excitatory Neurons/Progenitors,Excitatory Neurons/Progenitors,\
RGCs,RGCs,RGCs,RGCs,RGCs,RGCs,\
oRGCs,oRGCs,oRGCs,\
Pre-OPCs,Pre-OPCs,\
OPCs,OPCs,OPCs,OPCs,OPCs,OPCs,OPCs,OPCs,OPCs,OPCs,\
Pre-Astrocyte/Astrocytes,Pre-Astrocyte/Astrocytes,Pre-Astrocyte/Astrocytes,Pre-Astrocyte/Astrocytes,\
MGE Precursors,MGE Precursors,MGE Precursors,\
Cortical,Cortical,Cortical,Cortical,\
Striatal,Striatal,Striatal,Striatal,Striatal,\
SST+ (Cortex),SST+ (Cortex),SST+ (Cortex),SST+ (Cortex),\
LHX8+ (Cholinergic/Basal Telencephalon),LHX8+ (Cholinergic/Basal Telencephalon),LHX8+ (Cholinergic/Basal Telencephalon),LHX8+ (Cholinergic/Basal Telencephalon),\
CRABP1+,CRABP1+,CRABP1+,\
PV Precursors,PV Precursors,PV Precursors,PV Precursors,PV Precursors"

# ------------------------------------------------------------------------------
# Check that TARGET_GENES and TARGET_IDENTITIES have the same length
# ------------------------------------------------------------------------------
GENE_COUNT=$(echo "$TARGET_GENES" | tr ',' '\n' | wc -l)
IDENTITY_COUNT=$(echo "$TARGET_IDENTITIES" | tr ',' '\n' | wc -l)

echo "🧪 Checking counts..."
echo "  Number of genes:     $GENE_COUNT"
echo "  Number of identities: $IDENTITY_COUNT"

if [ "$GENE_COUNT" -ne "$IDENTITY_COUNT" ]; then
  echo "❌ Mismatch detected: TARGET_GENES and TARGET_IDENTITIES must have the same number of entries."
  exit 1
else
  echo "✅ Counts match. Proceeding with extraction."
fi

# Path to R script
SCRIPT_PATH="/home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/ExtractGenesByCluster.R"

# ------------------------------------------------------------------------------
# Run script
# ------------------------------------------------------------------------------
echo "📊 Extracting entries for multiple genes..."
Rscript "$SCRIPT_PATH"
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
  echo "❌ Extraction failed with exit code $EXIT_CODE"
  exit $EXIT_CODE
fi

echo "✅ Extraction complete"
