from mindscape.dagtoy.base_workflow import BaseWorkflow

class DataImportWorkflow(BaseWorkflow):
    def run(self):
        self.log_start()
        print(f"🚀 Running {self.__class__.__name__}")
        self.logger.info("Reached inside .run method of DataImportWorkflow")
        self.log_end()
