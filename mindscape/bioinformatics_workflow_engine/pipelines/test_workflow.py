from .base_workflow import BaseWorkflow

class TestWorkflow(BaseWorkflow):
    """
    TestWorkflow class that inherits from BaseWorkflow.
    Implements specific logic for the TestWorkflow pipeline.
    """

    def __init__(self, config):
        super().__init__(config)
        self.workflow_name = "TestWorkflow"
        self.setup_paths()

    def run(self):
        """
        Execute the TestWorkflow workflow.
        """
        if not self.config.get("force_rerun", False) and self.is_already_completed():
            print(f"✅ Skipping {self.workflow_name}; already completed.")
            return

        self.log_start()

        # TODO: Add your workflow steps here
        # Example: self.step_one()

        self.log_end()
        self.mark_completed()
