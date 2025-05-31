from pathlib import Path
import argparse
import yaml  # Ensure YAML is imported
from pipelines.base_workflow import BaseWorkflow
from pipelines.cell_ranger_workflow import CellRangerWorkflow
from pipelines.ventral_workflow import VentralWorkflow
from pipelines.qc_workflow import QCWorkflow
from utils.logger import setup_logger
from utils.config_merger import merge_configs

class WorkflowManager:
    def __init__(self, config_path, project_path):
        self.config_path = config_path  # Save the path to the configuration file
        self.config = self.load_config()
        self.logger = setup_logger("workflow_manager", "workflow_manager.log")
        self.workflows = []
        self.project_path = project_path  # Save the project path

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
                    # Pass the path to the configuration file instead of the loaded dictionary
                    self.workflows.append(workflow_class(config=self.config_path))
                else:
                    self.logger.warning(f"Workflow {workflow_name} not found.")

    def run_workflows(self):
        """Run all registered workflows using SLURM."""
        for workflow in self.workflows:
            workflow_name = workflow.__class__.__name__
            self.logger.info(f"Submitting workflow: {workflow_name} to SLURM")

            # Define the command to run the workflow with properly escaped quotes for embedding in a SLURM script
            command = (
                "python -c \""
                f"from pipelines.{workflow_name.lower()} import {workflow_name}; "
                f"{workflow_name}(\\\"{self.config_path}\\\").run()\""
            )

            # Submit the workflow as a SLURM job
            try:
                job_id = self.slurm_manager.submit_job(
                    command=command,
                    job_name=workflow_name,
                    pipeline_step="run"
                )
                self.logger.info(f"Workflow {workflow_name} submitted as SLURM job {job_id}")
            except RuntimeError as e:
                self.logger.error(f"Failed to submit workflow {workflow_name}: {e}")

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

    # Determine which configuration files to use
    if args.project_path:
        # If the user provides a `--project_path` argument, use it to locate the project-specific configuration file.
        project_path = Path(args.project_path)
        project_config_path = project_path / "config/config.yaml"
        if not project_config_path.exists():
            # Raise an error if the project-specific configuration file does not exist.
            raise FileNotFoundError(f"Configuration file not found at {project_config_path}")
        log_dir = project_path / "logs"  # Use the logs directory in the project path
    else:
        # If no `--project_path` is provided, fall back to default behavior.
        project_config_path = None
        log_dir = None  # Default to the current behavior

    # Define the path to the default configuration file
    default_config_path = Path(__file__).parent / "config/default_config.yaml"

    # Generate the merged configuration file
    if project_config_path:
        # If a project-specific configuration file exists, merge it with the default configuration.
        merged_config_path = merge_configs(default_config_path, project_config_path)
        print(f"Merged configuration saved to: {merged_config_path}")
    else:
        # If no project-specific configuration file exists, use the default configuration.
        merged_config_path = default_config_path

    # Initialize the WorkflowManager and run workflows
    logger = setup_logger("workflow_manager", "workflow_manager.log", log_dir=log_dir)
    workflow_manager = WorkflowManager(config_path=merged_config_path, project_path=args.project_path)
    workflow_manager.register_workflows()
    workflow_manager.run_workflows()