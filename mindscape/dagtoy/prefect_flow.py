from prefect import flow
from mindscape.dagtoy.create_project.new import create_new_project
from mindscape.dagtoy.runner import run_pipeline

@flow(name="MindScape Prefect DAG Test")
def run_mindscape_flow():
    project_path, config_path = create_new_project(
        base_dir="mindscape/dagtoy",
        project_name="TestDAGRun",
        experimenter="Manny"
    )
    run_pipeline(config_path)

if __name__ == "__main__":
    print("Loading MindScape 1.0.0...")
    run_mindscape_flow()