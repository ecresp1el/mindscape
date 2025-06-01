from .base_workflow import BaseWorkflow

class QCWorkflow(BaseWorkflow):
    """
    Quality Control Workflow for bioinformatics analysis.
    Inherits from BaseWorkflow to utilize shared logic.
    """

    def __init__(self, config):
        super().__init__(config)
        self.workflow_name = "QualityControlWorkflow"

    def run(self):
        """
        Execute the Quality Control workflow.
        This method will implement the specific steps for QC.
        """
        self.setup_paths()

        if not self.config.get("force_rerun", False) and self.is_already_completed():
            print(f"✅ Skipping {self.workflow_name}; already completed.")
            return

        self.log_start()

        # Implement QC steps here
        # Example: self.perform_qc_analysis()

        self.log_end()
        self.mark_completed()