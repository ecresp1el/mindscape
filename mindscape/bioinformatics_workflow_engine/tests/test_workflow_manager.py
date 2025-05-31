import unittest
import sys
from pathlib import Path

# Add the project root directory to sys.path
sys.path.insert(0, str(Path(__file__).parent.parent))

from pipelines.cell_ranger_workflow import CellRangerWorkflow
from pipelines.ventral_workflow import VentralWorkflow
from pipelines.qc_workflow import QCWorkflow

class TestWorkflowManager(unittest.TestCase):

    def setUp(self):
        # Provide an absolute path to the config file
        config_path = Path(__file__).parent.parent / 'config/default_config.yaml'
        self.cell_ranger_workflow = CellRangerWorkflow(config=config_path)
        self.ventral_workflow = VentralWorkflow(config=config_path)
        self.qc_workflow = QCWorkflow(config=config_path)

    def test_cell_ranger_workflow_initialization(self):
        self.assertIsInstance(self.cell_ranger_workflow, CellRangerWorkflow)

    def test_ventral_workflow_initialization(self):
        self.assertIsInstance(self.ventral_workflow, VentralWorkflow)

    def test_qc_workflow_initialization(self):
        self.assertIsInstance(self.qc_workflow, QCWorkflow)

    def test_workflow_run_method(self):
        # Test that the run method executes without errors
        try:
            self.cell_ranger_workflow.run()
            self.ventral_workflow.run()
            self.qc_workflow.run()
        except Exception as e:
            self.fail(f"Workflow run method raised an exception: {e}")

if __name__ == '__main__':
    unittest.main()