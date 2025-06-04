dag = {
    "AlignmentAndMoleculeCountingWorkflow": [],
    "CellFilteringWorkflow": ["AlignmentAndMoleculeCountingWorkflow"],
    "DoubletScoringWorkflow": ["CellFilteringWorkflow"],
    "CellSizeEstimationWorkflow": ["CellFilteringWorkflow"],
    "GeneVarianceAnalysisWorkflow": ["CellFilteringWorkflow"],
    "DimensionalityReductionWorkflow": ["GeneVarianceAnalysisWorkflow"],
    "ManifoldRepresentationWorkflow": ["DimensionalityReductionWorkflow"],
    "ClusteringAndDEWorkflow": ["ManifoldRepresentationWorkflow"],
    "TrajectoryInferenceWorkflow": ["DimensionalityReductionWorkflow"],
    "VelocityEstimationWorkflow": ["AlignmentAndMoleculeCountingWorkflow"],
    "CellTypeAnnotationWorkflow": ["ClusteringAndDEWorkflow"],
    "IntegrationWorkflow": ["GeneVarianceAnalysisWorkflow"],
    "MultiOmicsIntegrationWorkflow": ["IntegrationWorkflow"]
}
