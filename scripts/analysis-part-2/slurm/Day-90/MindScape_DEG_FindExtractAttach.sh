#!/bin/bash
#SBATCH --job-name=degs_find_extract_attach_day30
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=00:30:00
#SBATCH --output=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/degs_find_extract_attach_day30_%j.out
#SBATCH --error=/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/logs/degs_find_extract_attach_day30_%j.err

# ==============================================================================
# SLURM Job Script for MindScape - DEG Analysis + Extraction + Identity Append
# Combines Steps 3–5 of the MindScape pipeline.
# ==============================================================================

# -------------------- EMAIL NOTIFICATION ON FAILURE --------------------------
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oltij@umich.edu   # <-- change to your preferred address
# -----------------------------------------------------------------------------

echo "🚀 Starting MindScape DEG + Extract + Attach job"
date

# -------------------------- ENVIRONMENT SETUP ---------------------------------
module purge
source ~/.bashrc
conda activate mindscape-env

# -------------------------- PATH SETUP ----------------------------------------
export DEG_INPUT="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds"
export DEG_OUTPUT="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_clustered/degs_find_extract_attach_day30"
export NEW_INPUT_RDS="/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output_clustered_test_1/clustered_jeyoon_test_1.rds"
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

export MAX_CLUSTER_UNIQUENESS="13"

export TOP_N_PER_CLUSTER="5"

export OVERLAP_GENE_THRESHOLD=3

# -------------------------- CLUSTER IDENTITIES -------------------------------
# Define cluster identities as key-value pairs (to be parsed in R)
export CLUSTER_IDENTITIES="0=MGE Striatal/GP Fated,\
1=SST+, NPY +, Cortical Fated,\
2=CRABP1+/PV Precursors,\
3=PV precursors/Migrating cells/Cortical-fated,\
4=Pre-Astrocytes/Astrocytes 1,\
5=LHX8+ vMGE GABergic Striatal/GP fated 1,\
6=Stressed Cells,\
7=Stressed Cells,\
8=LHX8+ vMGE GABergic Striatal/GP fated 2,\
9=Pre-OPCs/OPCs,\
10=Pre-Astrocytes/Astrocytes 2,\
11=PV Precursors,\
12=Dividing cells"

# -------------------------- EXECUTION ----------------------------------------
Rscript /home/oltij/githubprojectfolder/mindscape/scripts/analysis-part-2/scripts/Day-90/MindScape_DEG_FindExtractAttach.R

status=$?
if [ $status -eq 0 ]; then
  echo "✅ MindScape DEG + Extract + Attach completed successfully"
else
  echo "❌ ERROR: Script exited with status $status"
fi

date
echo "🏁 Job finished."