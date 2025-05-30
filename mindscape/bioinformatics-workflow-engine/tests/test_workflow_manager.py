import unittest
from pipelines.base_workflow import BaseWorkflow
from pipelines.cell_ranger_workflow import CellRangerWorkflow
from pipelines.ventral_workflow import VentralWorkflow
from pipelines.qc_workflow import QCWorkflow
from utils.logger import setup_logger

class TestWorkflowManager(unittest.TestCase):

    def setUp(self):
        self.logger = setup_logger()
        self.cell_ranger_workflow = CellRangerWorkflow()
        self.ventral_workflow = VentralWorkflow()
        self.qc_workflow = QCWorkflow()

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
        self.cell_ranger_workflow.run()
        self.ventral_workflow.run()
        self.qc_workflow.run()
        # Add assertions to verify expected outcomes after running workflows

if __name__ == '__main__':
    unittest.main()