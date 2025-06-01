#!/usr/bin/env python
# NOTE: This script is located in tools/create_workflow_scaffold.py and is intended for developer use.
"""
🛠 MindScape Workflow Scaffold Generator

This script generates a new workflow module inside the `pipelines/` directory
using the MindScape modular workflow framework. It creates a new class that
inherits from `BaseWorkflow` and sets up boilerplate for configuration,
logging, and execution.

Usage (from project root):
    python tools/create_workflow_scaffold.py --name MyNewWorkflow

Options:
    --name WORKFLOW_NAME   Name of the workflow class and file (e.g., "CellRangerWorkflow")

Example:
    python tools/create_workflow_scaffold.py --name AnnotateClustersWorkflow

This will generate a new file at:
    mindscape/bioinformatics_workflow_engine/pipelines/annotate_clusters_workflow.py

And create a class:
    class AnnotateClustersWorkflow(BaseWorkflow):

Dependencies:
    - Must be run from the root of the MindScape project repo.
    - Assumes `BaseWorkflow` exists in `pipelines/base_workflow.py`.
"""
import argparse
import os
from pathlib import Path

TEMPLATE = '''from .base_workflow import BaseWorkflow

class {class_name}(BaseWorkflow):
    """
    {class_name} class that inherits from BaseWorkflow.
    Implements specific logic for the {class_name} pipeline.
    """

    def __init__(self, config):
        super().__init__(config)

    def run(self):
        """
        Execute the {class_name} workflow.
        """
        self.log_start()
        self.setup_paths()

        # TODO: Add your workflow steps here
        # Example: self.step_one()

        self.log_end()
'''

def to_snake_case(name):
    import re
    name = re.sub(r'(?<!^)(?=[A-Z])', '_', name).lower()
    return name

def main():
    parser = argparse.ArgumentParser(description="Generate a new MindScape workflow scaffold")
    parser.add_argument("--name", required=True, help="Name of the workflow class (e.g., MyNewWorkflow)")
    args = parser.parse_args()

    workflow_name = args.name
    class_name = workflow_name if workflow_name.endswith("Workflow") else workflow_name + "Workflow"
    filename = to_snake_case(class_name) + ".py"

    # Paths
    project_root = Path(__file__).resolve().parent.parent
    target_dir = project_root / "mindscape" / "bioinformatics_workflow_engine" / "pipelines"
    target_path = target_dir / filename

    if target_path.exists():
        print(f"❌ File already exists: {target_path}")
        return

    target_dir.mkdir(parents=True, exist_ok=True)
    with open(target_path, "w") as f:
        f.write(TEMPLATE.format(class_name=class_name))

    print(f"✅ Created workflow scaffold: {target_path}")
    print("📌 Next steps:")
    print(" - Implement your workflow logic inside the class")
    print(" - Add the workflow to your config.yaml under 'workflows'")
    print(" - Commit the new file to GitHub")
    print(" - This developer tool lives in: tools/create_workflow_scaffold.py")
    print(" - Run this script from the root of the MindScape repo")

if __name__ == "__main__":
    main()