import argparse
from pathlib import Path
from pipelines.base_workflow import BaseWorkflow
from utils.logger import setup_logger
import yaml

class WorkflowManager:
    def __init__(self, config_path):
        self.config_path = config_path
        self.config = self.load_config()
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
                    self.workflows.append(workflow_class(config=self.config))
                else:
                    self.logger.warning(f"Workflow {workflow_name} not found.")

    def run_workflows(self):
        """Run all registered workflows."""
        for workflow in self.workflows:
            self.logger.info(f"Starting workflow: {workflow.__class__.__name__}")
            workflow.run()
            self.logger.info(f"Completed workflow: {workflow.__class__.__name__}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Bioinformatics Workflows")
    parser.add_argument(
        "--project_path",
        type=str,
        required=False,
        help="Path to the project directory created by create_project",
    )
    args = parser.parse_args()

    # Determine which configuration file to use
    if args.project_path:
        project_path = Path(args.project_path)
        config_path = project_path / "config/config.yaml"
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found at {config_path}")
    else:
        # Fallback to default_config.yaml in the repository
        config_path = Path(__file__).parent / "config/default_config.yaml"

    # Load the configuration
    workflow_manager = WorkflowManager(config_path=config_path)
    workflow_manager.register_workflows()
    workflow_manager.run_workflows()