# mindscape/dagtoy/scripts/run_dag_test.py

from mindscape.dagtoy.create_project.new import create_new_project
import subprocess
from pathlib import Path
import sys

project_path, config_path = create_new_project(
    base_dir="mindscape/dagtoy", project_name="TestDAGRun", experimenter="Manny"
)

# Run the DAG using the generated config
command = [
    sys.executable,
    "mindscape/dagtoy/runner.py",
    "--config_path",
    str(config_path)
]

print("\n🚀 Launching DAG execution...")
subprocess.run(command, check=True)
print("✅ DAG run complete.")