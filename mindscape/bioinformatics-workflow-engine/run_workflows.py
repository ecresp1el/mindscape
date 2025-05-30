from pipelines.base_workflow import BaseWorkflow
from utils.logger import setup_logger
import yaml
from pathlib import Path

class WorkflowManager:
    def __init__(self, config_path: str):
        self.config_path = Path(config_path)
        self.workflows = []
        self.logger = setup_logger()

    def load_config(self):
        with open(self.config_path, 'r') as file:
            config = yaml.safe_load(file)
        return config

    def register_workflows(self, workflow_order):
        for workflow_name in workflow_order:
            module = __import__(f'pipelines.{workflow_name}', fromlist=[workflow_name])
            workflow_class = getattr(module, f'{workflow_name.capitalize()}Workflow')
            self.workflows.append(workflow_class())

    def run_workflows(self):
        for workflow in self.workflows:
            self.logger.info(f"Starting workflow: {workflow.__class__.__name__}")
            workflow.run()
            self.logger.info(f"Completed workflow: {workflow.__class__.__name__}")

if __name__ == "__main__":
    config = WorkflowManager('config/workflow_order.yaml')
    workflow_order = config.load_config().get('workflow_order', [])
    config.register_workflows(workflow_order)
    config.run_workflows()