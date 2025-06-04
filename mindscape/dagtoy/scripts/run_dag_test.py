# mindscape/dagtoy/scripts/run_dag_test.py

from mindscape.dagtoy.create_project.new import create_new_project
from mindscape.dagtoy.runner import run_pipeline
import sys 

# Create the test project
project_path, config_path = create_new_project(
    base_dir="mindscape/dagtoy", project_name="TestDAGRun", experimenter="Manny"
)

# Run the DAG directly
print("\n🚀 Launching DAG execution...")
run_pipeline(str(config_path))  # Important: pass config_path as string
print("✅ DAG run complete.")