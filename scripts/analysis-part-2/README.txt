MindScape — Day 90 Analysis Pipeline (mGEO SOSR Organoids)
=======================================================================

PURPOSE
-----------------------------------------------------------------------
This workflow performs integrated single-cell RNA-seq analysis for Day 90
mGEO organoid samples processed using 10x Genomics FLEX chemistry.

It takes normalized Seurat objects (from analysis-part-1) as inputs, performs integration,
clustering, and visualization of cell-type proportions via UMAP plots.

PIPELINE CONTEXT
-----------------------------------------------------------------------
This is the Day 90 branch of the MindScape pipeline, executed after:
  1. Cell Ranger preprocessing
  2. Sample-level normalization (Analysis Part 1)

STEPS AND SCRIPTS
-----------------------------------------------------------------------

1. **Integration (Part 1)**
   - Script: `mindscape_preprocess_integrated_jeyoon.R`
   - Purpose: Load normalized per-sample Seurat objects, merge and integrate
     them, regress out cell cycle effects, and save an integrated dataset
     for Day 90 analysis.

2. **Clustering (Part 2)**
   - Script: `mindscape_cluster_jeyoon.R`
   - Purpose: Run dimensionality reduction (PCA, UMAP) and clustering on
     the integrated Day 90 Seurat object. Outputs clustered RDS and basic plots.

3. **Marker Gene Identification (Part 3)**
   - Script: `MarkerGenes_AllClusters.R`
   - Purpose: Identify DEGs for each cluster vs. all others using Seurat’s
     `FindMarkers` function, producing verified per-cluster and combined CSVs.

4. **Gene Extraction (Part 4)**
   - Script: `ExtractGenesByCluster.R`
   - Purpose: Extract target genes of interest from the combined DEG table
     and organize them per cluster for downstream annotation or validation.

5. **UMAP + Cell-Type Proportion Visualization (Part 5)**
   - Script: `umap_proportions_figures.R`
   - Purpose: Generate verified UMAP and cluster proportion bar plots for
     the Day 90 dataset, combining both into final figure panels.

OUTPUTS
-----------------------------------------------------------------------
Each stage outputs verified RDS objects, CSVs, and plots within
the respective Day 90 output directories under:
   `/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/`

NOTES
-----------------------------------------------------------------------
- Run the scripts sequentially from Part 1 through Part 5.
- All scripts require the `mindscape-env` Conda environment
  (Seurat v5, dplyr, patchwork, ggplot2, Matrix, SingleCellExperiment).
- SLURM job scripts for each step are available under:
     `scripts/analysis-part-2/slurm/Day-90/`
- Verification and Annotation scripts can be utilized for ensuring
  indpendent runs match and can be annotated consistently using marker genes.

-----------------------------------------------------------------------
Contact:       oltij@umich.edu
Environment:   Great Lakes HPC — parent0 account
Pipeline:      MindScape (Day 90 — mGEO SOSR Organoids)
-----------------------------------------------------------------------
