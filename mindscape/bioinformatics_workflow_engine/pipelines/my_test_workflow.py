
from mindscape.bioinformatics_workflow_engine.baseworkflow import BaseWorkflow

class MyTestWorkflow(BaseWorkflow):
    def run(self):
        self.log_start()
        self.mark_in_progress()
        self.log_end()
        return True
