MindScape — mGEO Organoid Normalization Pipeline (FLEX 10x)
=======================================================================

PURPOSE
-----------------------------------------------------------------------
This step performs normalization and quality control on individual
10x Genomics FLEX scRNA-seq samples from Cell Ranger preprocessing.
Each sample is converted into a normalized Seurat object with QC metrics applied.

PIPELINE CONTEXT
-----------------------------------------------------------------------
This is the first step of the MindScape analysis pipeline, following
Cell Ranger output generation. The resulting normalized Seurat objects
are used as inputs for downstream integration, clustering, and DEG analysis.

Run order (Day 90):
1) mindscape_process_sample_normalized_jeyoon.R 
using run_ventral_sosrs_merged_jeyoon.sh