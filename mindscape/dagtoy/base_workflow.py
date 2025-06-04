from pathlib import Path

class BaseWorkflow:
    def __init__(self):
        self.name = self.__class__.__name__

    def is_complete(self):
        return Path(f"./test_files/{self.name}.done").exists()

    def mark_completed(self):
        Path("./test_files").mkdir(exist_ok=True)
        Path(f"./test_files/{self.name}.done").write_text("DONE")

    def run(self):
        print(f"🚀 Running {self.name}")
        self.mark_completed()
