# mindscape/dagtoy/base_workflow.py

from pathlib import Path
import logging
from datetime import datetime
import hashlib
import yaml


class BaseWorkflow:
    def __init__(self, config_path=None, logger=None, meta_hash=None):
        self.name = self.__class__.__name__

        if config_path is None:
            raise ValueError("❌ config_path must be provided to BaseWorkflow.")

        self.config_path = Path(config_path).resolve()
        self.logger = logger or logging.getLogger(f"workflow.{self.name}")
        self.logger.setLevel(logging.INFO)

        self.config_hash = self.compute_config_hash()
        self.meta_hash = meta_hash  # pipeline-level hash
        self.project_path = self.load_project_path_from_config()

        self.logfile = self.project_path / "logs" / f"{self.name}.log"
        self.logfile.parent.mkdir(exist_ok=True, parents=True)

        # Avoid adding duplicate handlers
        if not any(isinstance(h, logging.FileHandler) and getattr(h, "baseFilename", None) == str(self.logfile) for h in self.logger.handlers):
            file_handler = logging.FileHandler(self.logfile, mode='a')
            file_handler.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
            self.logger.addHandler(file_handler)

        if not any(isinstance(h, logging.StreamHandler) for h in self.logger.handlers):
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(logging.Formatter('%(message)s'))
            self.logger.addHandler(stream_handler)

        self.logger.debug(f"📄 Logging to file: {self.logfile}")

    def compute_config_hash(self):
        if not self.config_path.exists():
            return "NO_CONFIG"
        content = self.config_path.read_bytes()
        return hashlib.sha256(content).hexdigest()

    def is_complete(self):
        if not self.logfile.exists():
            return False
        lines = self.logfile.read_text().splitlines()
        for i in range(len(lines) - 1, -1, -1):
            if "Status: COMPLETED" in lines[i]:
                # look forward for config hash match
                context = "\n".join(lines[i:i+5])
                return self.config_hash in context
        return False

    def mark_completed(self):
        now = datetime.now().isoformat()
        self.logger.info(f"✅ Completed by: {self.__class__.__name__}")
        self.logger.info(f"Status: COMPLETED")
        self.logger.info(f"Timestamp: {now}")
        self.logger.info(f"Config hash: {self.config_hash}")
        if self.meta_hash:
            self.logger.info(f"Meta hash: {self.meta_hash}")

    # mark_in_progress is no longer needed; log_start handles status logging.

    def mark_failed(self, reason="Unspecified"):
        now = datetime.now().isoformat()
        self.logger.error(f"❌ Failed in: {self.__class__.__name__}")
        self.logger.error(f"Status: FAILED")
        self.logger.error(f"Timestamp: {now}")
        self.logger.error(f"Reason: {reason}")
        if self.meta_hash:
            self.logger.error(f"Meta hash: {self.meta_hash}")

    def get_status(self):
        # Optionally update this method to check log tail for statuses
        if not self.logfile.exists():
            return "not_started"
        tail = self.logfile.read_text().splitlines()[-10:]
        joined = "\n".join(tail)
        if "Status: COMPLETED" in joined and self.config_hash in joined:
            return "completed"
        if "Status: FAILED" in joined:
            return "failed"
        if "Status: IN PROGRESS" in joined:
            return "in_progress"
        return "not_started"

    # get_marker is no longer needed; marker files are not used.

    def log_start(self):
        now = datetime.now().isoformat()
        self.logger.info(f"🟢 Starting workflow: {self.name}")
        self.logger.info(f"Status: IN PROGRESS")
        self.logger.info(f"Timestamp: {now}")
        self.logger.info(f"Config hash: {self.config_hash}")
        if self.meta_hash:
            self.logger.info(f"Meta hash: {self.meta_hash}")

    def log_end(self):
        self.logger.info(f"✅ Completed workflow: {self.name}")
        self.mark_completed()

    def run(self):
        raise NotImplementedError("Subclasses must implement the run() method.")


    def load_project_path_from_config(self):
        if not self.config_path.exists():
            raise FileNotFoundError(f"❌ config_path does not exist: {self.config_path}")
        with open(self.config_path, "r") as f:
            config = yaml.safe_load(f)
        if "project_path" not in config:
            raise KeyError("❌ 'project_path' key is missing from config file.")
        return Path(config["project_path"]).resolve()