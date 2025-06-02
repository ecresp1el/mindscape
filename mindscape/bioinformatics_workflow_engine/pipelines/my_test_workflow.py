
from .base_workflow import BaseWorkflow

class MyTestWorkflow(BaseWorkflow):
    def run(self):
        self.log_start()
        self.mark_in_progress()
        self.log_end()
        return True

# Simulated failure workflow for testing mark_failed()
class FailingWorkflow(BaseWorkflow):
    def run(self):
        self.log_start()
        self.mark_in_progress()
        raise RuntimeError("Simulated workflow failure for testing mark_failed()")
