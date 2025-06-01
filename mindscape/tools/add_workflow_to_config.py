# add_workflow_to_config.py
"""
📘 Add Workflow to MindScape config.yaml

This developer tool adds a new workflow to an existing MindScape project configuration.
It ensures:
  - The workflow is only added if it doesn't already exist.
  - The workflow is enabled by default (unless --disabled is used).
  - The new workflow name will also be dynamically importable if added to config.
  - Warns if the workflow file does not exist in the expected pipeline folder.

Usage:
------
Run from the root of the MindScape GitHub repository:

    python mindscape/tools/add_workflow_to_config.py \
        --project_path /nfs/turbo/umms-parent/myproject-MannyCrespo-2025-06-01 \
        --workflow_name MyCustomWorkflow

Optional flags:
  --disabled       Add the workflow with "enabled: false"

This script does **not** create the workflow Python file (see `create_workflow_scaffold.py` for that).
It only modifies the config.yaml and does not register the workflow class itself.
"""

import argparse
import os
from pathlib import Path
from ruamel.yaml import YAML


def add_workflow_to_config(config_path: Path, workflow_name: str, enabled: bool = True):
    yaml = YAML()
    yaml.preserve_quotes = True

    if not config_path.exists():
        raise FileNotFoundError(f"❌ Config file not found: {config_path}")

    with config_path.open("r") as f:
        config = yaml.load(f)

    if "workflows" not in config:
        config["workflows"] = []

    existing_names = {wf.get("name") for wf in config["workflows"]}
    pipeline_file = Path(__file__).parent.parent / "bioinformatics_workflow_engine" / "pipelines" / f"{workflow_name.lower()}.py"
    if not pipeline_file.exists():
        print(f"⚠️ Warning: No Python file found for workflow '{workflow_name}' at {pipeline_file}")
    if workflow_name in existing_names:
        print(f"⚠️ Workflow '{workflow_name}' already exists in config. Skipping.")
        return

    config["workflows"].append({"name": workflow_name, "enabled": enabled})

    with config_path.open("w") as f:
        yaml.dump(config, f)

    print(f"✅ Added workflow '{workflow_name}' to config (enabled={enabled})")


if __name__ == "__main__":
    # This only modifies the config.yaml and does not register the workflow class itself.
    parser = argparse.ArgumentParser(description="Add a workflow to an existing MindScape project config.")
    parser.add_argument("--project_path", type=str, required=True, help="Path to the MindScape project directory")
    parser.add_argument("--workflow_name", type=str, required=True, help="Name of the workflow class (must match .py filename)")
    parser.add_argument("--disabled", action="store_true", help="Add the workflow as disabled")
    args = parser.parse_args()

    project_config_path = Path(args.project_path) / "config" / "config.yaml"
    add_workflow_to_config(project_config_path, args.workflow_name, enabled=not args.disabled)