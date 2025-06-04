#!/bin/bash

# Always run from this script's directory
cd "$(dirname "$0")" || { echo "❌ Failed to enter script directory."; exit 1; }

# Setup folders
mkdir -p workflows test_files
touch __init__.py workflows/__init__.py

# Helper function to convert CamelCase to snake_case
camel_to_snake() {
    echo "$1" | sed -E 's/([a-z])([A-Z])/\1_\2/g' | sed -E 's/([A-Z]+)([A-Z][a-z])/\1_\2/g' | tr '[:upper:]' '[:lower:]'
}

# Generate all workflow files
workflow_classes=(
    DataImportWorkflow
    AlignmentAndMoleculeCountingWorkflow
    CellFilteringWorkflow
    DoubletScoringWorkflow
    CellSizeEstimationWorkflow
    GeneVarianceAnalysisWorkflow
    DimensionalityReductionWorkflow
    ManifoldRepresentationWorkflow
    ClusteringAndDEWorkflow
    TrajectoryInferenceWorkflow
    VelocityEstimationWorkflow
    CellTypeAnnotationWorkflow
    IntegrationWorkflow
    MultiOmicsIntegrationWorkflow
)

for class_name in "${workflow_classes[@]}"; do
    filename=$(camel_to_snake "${class_name%Workflow}").py
    cat > "workflows/$filename" <<EOF
from mindscape.dagtoy.base_workflow import BaseWorkflow

class ${class_name}(BaseWorkflow):
    def run(self):
        print(f"🚀 Running {self.__class__.__name__}")
        self.mark_completed()
EOF
done

echo "✅ DAGToy scaffold created inside: $(pwd)"