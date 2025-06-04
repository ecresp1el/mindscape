# mindscape/dagtoy/base_workflow.py

from pathlib import Path
import logging

class BaseWorkflow:
    def __init__(self, config_path="dagtoy/test_config.yaml", logger=None):
        self.name = self.__class__.__name__
        self.config_path = Path(config_path)
        self.project_path = Path("dagtoy/test_files")
        self.logger = logger or logging.getLogger(f"workflow.{self.name}")
        self.logger.setLevel(logging.INFO)
        self.logfile = self.project_path / f"{self.name}.log"
        self.logfile.parent.mkdir(exist_ok=True, parents=True)

    def is_complete(self):
        return self.get_marker(".done").exists()

    def mark_completed(self):
        self.get_marker(".done").write_text("COMPLETED\n")

    def mark_in_progress(self):
        self.get_marker(".in_progress").write_text("IN PROGRESS\n")

    def mark_failed(self, reason="Unspecified"):
        self.get_marker(".failed").write_text(f"FAILED: {reason}\n")

    def get_status(self):
        if self.get_marker(".done").exists():
            return "completed"
        if self.get_marker(".in_progress").exists():
            return "in_progress"
        if self.get_marker(".failed").exists():
            return "failed"
        return "not_started"

    def get_marker(self, suffix):
        return self.project_path / f"{self.name}{suffix}"

    def log_start(self):
        self.logger.info(f"🟢 Starting workflow: {self.name}")
        self.mark_in_progress()

    def log_end(self):
        self.logger.info(f"✅ Completed workflow: {self.name}")
        self.mark_completed()

    def run(self):
        raise NotImplementedError("Subclasses must implement the run() method.")

