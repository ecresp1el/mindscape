from mindscape.dagtoy.base_workflow import BaseWorkflow

class CellFilteringWorkflow(BaseWorkflow):
    def run(self):
        print(f"🚀 Running {self.__class__.__name__}")
        self.mark_completed()
