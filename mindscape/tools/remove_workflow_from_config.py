"""
Remove Workflow from MindScape config.yaml

This developer tool removes a workflow entry from an existing MindScape project configuration.
It only affects the `config.yaml` file.

Usage:
------
Run from the root of the MindScape GitHub repository:

    python mindscape/tools/remove_workflow_from_config.py \
        --project_path /nfs/turbo/umms-parent/myproject-MannyCrespo-2025-06-01 \
        --workflow_name MyCustomWorkflow

This script does **not** delete the Python file or unregister the class.
It simply edits the config.
"""

import argparse
from pathlib import Path
from ruamel.yaml import YAML

def remove_workflow_from_config(config_path: Path, workflow_name: str):
    yaml = YAML()
    yaml.preserve_quotes = True

    if not config_path.exists():
        raise FileNotFoundError(f"❌ Config file not found: {config_path}")

    with config_path.open("r") as f:
        config = yaml.load(f)

    if "workflows" not in config or not isinstance(config["workflows"], list):
        print("⚠️ No workflows section found in config.")
        return

    original_len = len(config["workflows"])
    config["workflows"] = [wf for wf in config["workflows"] if wf.get("name") != workflow_name]

    if len(config["workflows"]) == original_len:
        print(f"⚠️ Workflow '{workflow_name}' not found in config.")
        return

    with config_path.open("w") as f:
        yaml.dump(config, f)

    print(f"✅ Removed workflow '{workflow_name}' from config.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove a workflow from a MindScape config.yaml")
    parser.add_argument("--project_path", type=str, required=True, help="Path to the MindScape project directory")
    parser.add_argument("--workflow_name", type=str, required=True, help="Name of the workflow to remove")
    args = parser.parse_args()

    config_path = Path(args.project_path) / "config" / "config.yaml"
    remove_workflow_from_config(config_path, args.workflow_name)

