import os
import re
from pathlib import Path
import argparse
import yaml  # Ensure YAML is imported
from pipelines.base_workflow import BaseWorkflow
from pipelines.cell_ranger_workflow import CellRangerWorkflow
from pipelines.ventral_workflow import VentralWorkflow
from pipelines.qc_workflow import QCWorkflow
from utils.logger import setup_logger
from utils.validation import warn_if_missing_from_config


class WorkflowManager:
    def __init__(self, config_path, project_path):
        self.config_path = config_path  # Save the path to the configuration file
        self.project_path = project_path  # Save the project path
        self.config = self.load_config()
        # Warn about workflows that exist as files but aren't declared in config.yaml
        pipeline_dir = Path(__file__).parent / "pipelines"
        configured_names = {wf["name"] for wf in self.config.get("workflows", [])}
        warn_if_missing_from_config(pipeline_dir, configured_names)
        self.logger = setup_logger("workflow_manager", "workflow_manager.log")
        self.workflows = []

    def load_config(self):
        """Load the YAML configuration file."""
        with open(self.config_path, 'r') as file:
            return yaml.safe_load(file)

    def register_workflows(self):
        """Register workflows based on the configuration."""
        workflow_order = self.config.get("workflows", [])
        for workflow in workflow_order:
            workflow_name = workflow.get("name")
            if workflow.get("enabled", False):
                workflow_class = globals().get(workflow_name)
                if workflow_class:
                    # Check if the workflow has overridden the BaseWorkflow.run() method
                    # This prevents registration of placeholder or scaffolded workflows
                    #
                    # Explanation:
                    # - If a workflow class does not implement its own .run() method, it inherits BaseWorkflow.run.
                    # - Registering such a workflow is likely a mistake (scaffolded but not implemented).
                    # - This check skips registration and warns the developer.
                    if workflow_class.run == BaseWorkflow.run:
                        print(f"⚠️ Skipping '{workflow_name}': .run() is not implemented.")
                        continue
                    self.workflows.append(workflow_class(config=self.config_path))
                else:
                    self.logger.warning(f"Workflow {workflow_name} not found.")

    def run_workflows(self):
        """Run all registered workflows with per-step completion checks."""
        force_rerun = self.config.get("force_rerun", False)

        for workflow in self.workflows:
            workflow_name = workflow.__class__.__name__
            self.logger.info(f"Running workflow: {workflow_name}")

            # Only apply per-workflow skipping if force_rerun is False
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
    # This line checks if the script is being run directly (as the main program) or being imported as a module.
    #
    # Explanation:
    # - When a Python script is executed, the special variable `__name__` is automatically set by Python.
    # - If the script is being run directly (e.g., `python run_workflows.py`), `__name__` is set to `"__main__"`.
    # - If the script is being imported as a module into another script, `__name__` is set to the name of the module
    #   (e.g., `"run_workflows"`).
    #
    # Why is this check important?
    # - This check ensures that the code inside this block is only executed when the script is run directly.
    # - If the script is imported as a module, the code inside this block will NOT run, preventing unintended behavior.
    #
    # In this case, the block initializes the argument parser, determines the configuration file paths, and
    # runs the workflows. This is the main entry point for the script.

    parser = argparse.ArgumentParser(description="Run Bioinformatics Workflows")
    parser.add_argument(
        "--project_path",
        type=str,
        required=False,
        help="Path to the project directory created by create_project",
    )
    args = parser.parse_args()

    parser = argparse.ArgumentParser(description="Run Bioinformatics Workflows")
    parser.add_argument(
        "--project_path",
        type=str,
        required=True,
        help="Path to the project directory created by create_project",
    )
    args = parser.parse_args()

    project_path = Path(args.project_path)
    project_config_path = project_path / "config/config.yaml"
    if not project_config_path.exists():
        raise FileNotFoundError(f"Configuration file not found at {project_config_path}")
    log_dir = project_path / "logs"

    logger = setup_logger("workflow_manager", "workflow_manager.log", log_dir=log_dir)
    workflow_manager = WorkflowManager(config_path=project_config_path, project_path=args.project_path)
    workflow_manager.register_workflows()
    workflow_manager.run_workflows()