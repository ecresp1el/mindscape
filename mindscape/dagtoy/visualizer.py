# mindscape/dagtoy/visualizer.py
from graphviz import Digraph
from mindscape.dagtoy.dag_config import dag
import os

# Define node statuses manually (used for color-coding)
NODE_STATUS = {
    "DataImportWorkflow": "real",
    "AlignmentAndMoleculeCountingWorkflow": "real",
    "CellFilteringWorkflow": "real",
    "DoubletScoringWorkflow": "placeholder",
    "CellSizeEstimationWorkflow": "placeholder",
    "GeneVarianceAnalysisWorkflow": "real",
    "DimensionalityReductionWorkflow": "real",
    "ManifoldRepresentationWorkflow": "real",
    "ClusteringAndDEWorkflow": "real",
    "TrajectoryInferenceWorkflow": "placeholder",
    "VelocityEstimationWorkflow": "placeholder",
    "CellTypeAnnotationWorkflow": "placeholder",
    "IntegrationWorkflow": "real",
    "MultiOmicsIntegrationWorkflow": "placeholder"
}

COLOR_MAP = {
    "real": "lightgreen",
    "placeholder": "lightgray",
    "test": "lightblue"
}


def visualize_dag(output_path="mindscape/dagtoy/outputs/dag_structure", filetype="png"):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    dot = Digraph(comment="MindScape Toy DAG", format=filetype)
    dot.attr(rankdir="TB", size="8,10")

    for step in dag:
        label = step.replace("Workflow", "")
        status = NODE_STATUS.get(step, "placeholder")
        color = COLOR_MAP.get(status, "white")
        dot.node(step, label=label, style="filled", fillcolor=color, shape="box")

    for step, deps in dag.items():
        for dep in deps:
            dot.edge(dep, step)

    dot.render(output_path, cleanup=True)
    print(f"✅ DAG visual written to {output_path}.{filetype}")

if __name__ == "__main__":
    visualize_dag(output_path="mindscape/dagtoy/outputs/dag_structure", filetype="png")