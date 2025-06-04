#!/bin/bash

# Navigate to dagtoy directory
cd dagtoy || { echo "❌ dagtoy directory not found."; exit 1; }

# Create folders
mkdir -p workflows test_files

# Init files
touch __init__.py
touch workflows/__init__.py

# Base workflow stub
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

# DAG configuration
cat > dag_config.py <<EOF
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
EOF

# DAG runner
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

def run_pipeline():
    order = topological_sort(dag)
    for name in order:
        module_name = f"dagtoy.workflows.{name[:-9].lower()}"
        class_name = name
        module = importlib.import_module(module_name)
        klass = getattr(module, class_name)
        instance = klass()
        if not instance.is_complete():
            instance.run()
        else:
            print(f"✅ Skipping {class_name} (already complete)")

if __name__ == "__main__":
    run_pipeline()
EOF

# 13 workflow stubs
for name in \
    AlignmentAndMoleculeCounting \
    CellFiltering \
    DoubletScoring \
    CellSizeEstimation \
    GeneVarianceAnalysis \
    DimensionalityReduction \
    ManifoldRepresentation \
    ClusteringAndDE \
    TrajectoryInference \
    VelocityEstimation \
    CellTypeAnnotation \
    Integration \
    MultiOmicsIntegration
do
cat > workflows/${name,,}.py <<EOF
from dagtoy.base_workflow import BaseWorkflow

class ${name}Workflow(BaseWorkflow):
    def run(self):
        print(f"🚀 Running {self.__class__.__name__}")
        self.mark_completed()
EOF
done

echo "✅ DAGToy scaffold created!"