# mindscape/dagtoy/prefect_flow.py
from prefect import flow
from mindscape.dagtoy.create_project.new import create_new_project
from mindscape.dagtoy.workflows import *



@flow(name="MindScape Prefect DAG", log_prints=True)
def run_mindscape_flow():
    # Create project structure and get config paths
    project_path, _ = create_new_project(
        base_dir="mindscape/dagtoy",
        project_name="TestDAGRun",
        experimenter="Manny"
    )

    # Instantiate and run each workflow step
    DataImportWorkflow(project_path).run()
    AlignmentAndMoleculeCountingWorkflow(project_path).run()
    CellFilteringWorkflow(project_path).run()
    DoubletScoringWorkflow(project_path).run()
    CellSizeEstimationWorkflow(project_path).run()
    GeneVarianceAnalysisWorkflow(project_path).run()
    DimensionalityReductionWorkflow(project_path).run()
    ManifoldRepresentationWorkflow(project_path).run()
    ClusteringAndDEWorkflow(project_path).run()
    TrajectoryInferenceWorkflow(project_path).run()
    VelocityEstimationWorkflow(project_path).run()
    CellTypeAnnotationWorkflow(project_path).run()
    IntegrationWorkflow(project_path).run()
    MultiOmicsIntegrationWorkflow(project_path).run()


