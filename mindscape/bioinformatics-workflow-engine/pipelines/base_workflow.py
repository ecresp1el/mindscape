from __future__ import annotations

import os
import logging
from pathlib import Path
import yaml

class BaseWorkflow:
    def __init__(self, config_path: str | Path):
        self.config_path = Path(config_path)
        self.config = self.load_config()
        self.setup_logging()

    def load_config(self) -> dict:
        """Load the YAML configuration file."""
        if not self.config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.config_path}")
        with open(self.config_path, 'r') as file:
            return yaml.safe_load(file)

    def setup_logging(self):
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(self.__class__.__name__)

    def set_up_paths(self, project_path):
        """Set up project paths."""
        self.project_path = Path(project_path)
        if not self.project_path.exists():
            self.project_path.mkdir(parents=True, exist_ok=True)

    def run(self):
        """Run the workflow. Must be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement the run method.")