# mindscape/dagtoy/base_workflow.py

from pathlib import Path
import yaml
from datetime import datetime
import hashlib

class BaseWorkflow:
    def __init__(self, config_path: str):
        self.name = self.__class__.__name__
        self.config_path = Path(config_path).resolve()
        self.project_path = self._load_project_path()
        self.config_hash = self._compute_config_hash()

    def _compute_config_hash(self):
        content = self.config_path.read_bytes()
        return hashlib.sha256(content).hexdigest()

    def _load_project_path(self):
        if not self.config_path.exists():
            raise FileNotFoundError(f"❌ config_path does not exist: {self.config_path}")
        with open(self.config_path, "r") as f:
            config = yaml.safe_load(f)
        if "project_path" not in config:
            raise KeyError("❌ 'project_path' is missing from config file.")
        return Path(config["project_path"]).resolve()

    def is_complete(self):
        # Placeholder for future hash-based or file-based completion check
        return False

    def run(self):
        raise NotImplementedError("Subclasses must implement the run() method.")