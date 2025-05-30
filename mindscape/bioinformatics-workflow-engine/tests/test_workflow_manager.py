import unittest
from pipelines.cell_ranger_workflow import CellRangerWorkflow
from pipelines.ventral_workflow import VentralWorkflow
from pipelines.qc_workflow import QCWorkflow
from utils.logger import setup_logger

class TestWorkflowManager(unittest.TestCase):

    def setUp(self):
        self.logger = setup_logger("test_workflow_manager", "test_workflow_manager.log")
        test_config_path = 'config/default_config.yaml'
        self.cell_ranger_workflow = CellRangerWorkflow(config=test_config_path)
        self.ventral_workflow = VentralWorkflow(config=test_config_path)
        self.qc_workflow = QCWorkflow(config=test_config_path)

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