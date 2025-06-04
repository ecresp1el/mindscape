#!/bin/bash

# Always run from this script's directory
cd "$(dirname "$0")" || { echo "❌ Failed to enter script directory."; exit 1; }

# Setup folders
mkdir -p workflows test_files
touch __init__.py workflows/__init__.py

# base_workflow.py
cat > base_workflow.py <<EOF
# mindscape/dagtoy/base_workflow.py

from pathlib import Path
import logging

class BaseWorkflow:
    def __init__(self, config_path="dagtoy/test_config.yaml", logger=None):
        self.name = self.__class__.__name__
        self.config_path = Path(config_path)
        self.project_path = Path("dagtoy/test_files")
        self.logger = logger or logging.getLogger(f"workflow.{self.name}")
        self.logger.setLevel(logging.INFO)
        self.logfile = self.project_path / f"{self.name}.log"
        self.logfile.parent.mkdir(exist_ok=True, parents=True)

    def is_complete(self):
        return self.get_marker(".done").exists()

    def mark_completed(self):
        self.get_marker(".done").write_text("COMPLETED\n")

    def mark_in_progress(self):
        self.get_marker(".in_progress").write_text("IN PROGRESS\n")

    def mark_failed(self, reason="Unspecified"):
        self.get_marker(".failed").write_text(f"FAILED: {reason}\n")

    def get_status(self):
        if self.get_marker(".done").exists():
            return "completed"
        if self.get_marker(".in_progress").exists():
            return "in_progress"
        if self.get_marker(".failed").exists():
            return "failed"
        return "not_started"

    def get_marker(self, suffix):
        return self.project_path / f"{self.name}{suffix}"

    def log_start(self):
        self.logger.info(f"🟢 Starting workflow: {self.name}")
        self.mark_in_progress()

    def log_end(self):
        self.logger.info(f"✅ Completed workflow: {self.name}")
        self.mark_completed()

    def run(self):
        raise NotImplementedError("Subclasses must implement the run() method.")

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