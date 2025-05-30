from .base_workflow import BaseWorkflow

class CellRangerWorkflow(BaseWorkflow):
    def __init__(self, config):
        super().__init__(config)
        self.workflow_name = "Cell Ranger Workflow"

    def run(self):
        """Execute the Cell Ranger workflow steps."""
        self.setup_paths()  # This now works as expected
        self.step1()
        self.step2()
        self.step3()

    def step1(self):
        """Placeholder for the first step of the Cell Ranger workflow."""
        pass

    def step2(self):
        """Placeholder for the second step of the Cell Ranger workflow."""
        pass

    def step3(self):
        """Placeholder for the third step of the Cell Ranger workflow."""
        pass