#!/bin/bash

# Always run from this script's directory
cd "$(dirname "$0")" || { echo "❌ Failed to enter script directory."; exit 1; }

# Setup folders
mkdir -p workflows test_files
touch __init__.py workflows/__init__.py

# base_workflow.py
cat > base_workflow.py <<EOF
from pathlib import Path

class BaseWorkflow:
    def __init__(self):
        self.name = self.__class__.__name__

    def is_complete(self):
        return Path(f"./test_files/{self.name}.done").exists()

    def mark_completed(self):
        Path("./test_files").mkdir(exist_ok=True)
        Path(f"./test_files/{self.name}.done").write_text("DONE")

    def run(self):
        print(f"🚀 Running {self.name}")
        self.mark_completed()
EOF

# dag_config.py
cat > dag_config.py <<EOF
dag = {
    
    "DataImportWorkflow": [],
    "AlignmentAndMoleculeCountingWorkflow": ["DataImportWorkflow"],
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
EOF

# runner.py with removesuffix + snake_case conversion
cat > runner.py <<EOF
import importlib
from dag_config import dag

def topological_sort(graph):
    visited, order = set(), []
    def visit(node):
        if node not in visited:
            for dep in graph.get(node, []):
                visit(dep)
            visited.add(node)
            order.append(node)
    for node in graph:
        visit(node)
    return order

def snake_case(name):
    import re
    s1 = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\\1_\\2", s1).lower()

def run_pipeline():
    order = topological_sort(dag)
    for name in order:
        if not name.endswith("Workflow"):
            raise ValueError(f"Workflow name '{name}' must end in 'Workflow'")

        base = name.removesuffix("Workflow")
        snake = snake_case(base)
        module = importlib.import_module(f"mindscape.dagtoy.workflows.{snake}")
        klass = getattr(module, name)
        instance = klass()

        if not instance.is_complete():
            instance.run()
        else:
            print(f"✅ Skipping {name} (already complete)")

if __name__ == "__main__":
    run_pipeline()
EOF

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
echo "▶️  Running test pipeline..."

PYTHONPATH=../.. python runner.py