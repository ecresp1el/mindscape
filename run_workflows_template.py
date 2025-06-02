"""
🚀 Blank WorkflowManager Template Generator

This script generates a blank `run_workflows.py` file that:
- Includes only BaseWorkflow
- Does NOT pre-wire any existing subclasses
- Is ready for dynamic registration via config

You can then manually add new workflow imports or dynamically scan them.
"""

import os
import re
from pathlib import Path
import argparse
import yaml

from pipelines.base_workflow import BaseWorkflow
from utils.logger import setup_logger
from utils.validation import warn_if_missing_from_config



class WorkflowManager:
    def __init__(self, config_path, project_path):
        self.config_path = config_path
        self.project_path = project_path
        self.config = self.load_config()
        pipeline_dir = Path(__file__).parent / "pipelines"
        configured_names = {wf["name"] for wf in self.config.get("workflows", [])}
        warn_if_missing_from_config(pipeline_dir, configured_names)
        self.logger = setup_logger("workflow_manager", "workflow_manager.log")
        self.workflows = []

    def load_config(self):
        with open(self.config_path, 'r') as file:
            return yaml.safe_load(file)

    def register_workflows(self):
        workflow_order = self.config.get("workflows", [])
        for workflow in workflow_order:
            workflow_name = workflow.get("name")
            if workflow.get("enabled", False):
                workflow_class = globals().get(workflow_name)
                if workflow_class and hasattr(workflow_class, "run") and workflow_class.run != BaseWorkflow.run:
                    self.workflows.append(workflow_class(config=self.config_path))
                else:
                    self.logger.warning(f"Workflow {workflow_name} not found or not implemented.")

    def run_workflows(self):
        force_rerun = self.config.get("force_rerun", False)
        for workflow in self.workflows:
            workflow_name = workflow.__class__.__name__
            self.logger.info(f"Running workflow: {workflow_name}")
            if not force_rerun:
                if hasattr(workflow, "is_already_completed") and workflow.is_already_completed():
                    self.logger.info(f"✅ Skipping {workflow_name}: already completed.")
                    continue
            try:
                workflow.run()
                if hasattr(workflow, "mark_completed"):
                    workflow.mark_completed()
                self.logger.info(f"✅ Workflow {workflow_name} completed successfully.")
            except RuntimeError as e:
                self.logger.error(f"❌ Failed to run workflow {workflow_name}: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Bioinformatics Workflows")
    parser.add_argument("--project_path", type=str, required=True, help="Path to the MindScape project")
    args = parser.parse_args()

    project_path = Path(args.project_path)
    config_path = project_path / "config/config.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found at {config_path}")

    logger = setup_logger("workflow_manager", "workflow_manager.log", log_dir=project_path / "logs")
    workflow_manager = WorkflowManager(config_path=config_path, project_path=project_path)
    workflow_manager.register_workflows()
    workflow_manager.run_workflows()
