# mindscape/dagtoy/base_workflow.py

from pathlib import Path
import logging
from datetime import datetime
import hashlib

class BaseWorkflow:
    def __init__(self, config_path="dagtoy/test_config.yaml", logger=None):
        self.name = self.__class__.__name__
        self.config_path = Path(config_path)
        self.project_path = Path("dagtoy/test_files")
        self.logger = logger or logging.getLogger(f"workflow.{self.name}")
        self.logger.setLevel(logging.INFO)
        self.logfile = self.project_path / f"{self.name}.log"
        self.logfile.parent.mkdir(exist_ok=True, parents=True)
        self.config_hash = self.compute_config_hash()

    def compute_config_hash(self):
        if not self.config_path.exists():
            return "NO_CONFIG"
        content = self.config_path.read_bytes()
        return hashlib.sha256(content).hexdigest()

    def is_complete(self):
        marker = self.get_marker(".done")
        if not marker.exists():
            return False
        if self.config_hash not in marker.read_text():
            self.logger.info(f"⚠️ Config changed since last run of {self.name}")
            return False
        return True

    def mark_completed(self):
        now = datetime.now().isoformat()
        text = f"COMPLETED at {now}\nCONFIG_HASH: {self.config_hash}\n"
        self.get_marker(".done").write_text(text)

    def mark_in_progress(self):
        now = datetime.now().isoformat()
        self.get_marker(".in_progress").write_text(f"IN PROGRESS at {now}\n")

    def mark_failed(self, reason="Unspecified"):
        now = datetime.now().isoformat()
        self.get_marker(".failed").write_text(f"FAILED at {now}\nREASON: {reason}\n")

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

