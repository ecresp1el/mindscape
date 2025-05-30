import unittest
from pipelines.base_workflow import BaseWorkflow
from pipelines.cell_ranger_workflow import CellRangerWorkflow
from pipelines.ventral_workflow import VentralWorkflow
from pipelines.qc_workflow import QCWorkflow
from utils.logger import setup_logger

class TestWorkflowManager(unittest.TestCase):

    def setUp(self):
        # Provide valid arguments for setup_logger
        self.logger = setup_logger("test_workflow_manager", "test_workflow_manager.log")
        
        # Provide a valid config path for testing
        test_config_path = 'config/default_config.yaml'
        
        # Initialize workflows with the required config argument
        self.cell_ranger_workflow = CellRangerWorkflow(config=test_config_path)
        self.ventral_workflow = VentralWorkflow(config=test_config_path)
        self.qc_workflow = QCWorkflow(config=test_config_path)

    def test_cell_ranger_workflow_initialization(self):
        self.assertIsInstance(self.cell_ranger_workflow, CellRangerWorkflow)
        self.assertIsInstance(self.cell_ranger_workflow, BaseWorkflow)

    def test_ventral_workflow_initialization(self):
        self.assertIsInstance(self.ventral_workflow, VentralWorkflow)
        self.assertIsInstance(self.ventral_workflow, BaseWorkflow)

    def test_qc_workflow_initialization(self):
        self.assertIsInstance(self.qc_workflow, QCWorkflow)
        self.assertIsInstance(self.qc_workflow, BaseWorkflow)

    def test_workflow_run_method(self):
        # Ensure the run method works for each workflow
        with self.assertRaises(NotImplementedError):
            self.cell_ranger_workflow.run()
        with self.assertRaises(NotImplementedError):
            self.ventral_workflow.run()
        with self.assertRaises(NotImplementedError):
            self.qc_workflow.run()

if __name__ == '__main__':
    unittest.main()