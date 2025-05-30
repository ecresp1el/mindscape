from __future__ import annotations

import os
import logging
from pathlib import Path
import yaml

class BaseWorkflow:
    def __init__(self, config_path):
        # Resolve the config path to an absolute path
        self.config_path = Path(config_path).resolve()
        self.config = self.load_config()
        self.workflow_name = self.__class__.__name__

    def load_config(self):
        """Load the YAML configuration file."""
        with open(self.config_path, 'r') as file:
            return yaml.safe_load(file)

    def setup_paths(self):
        """Set up project paths."""
        self.project_path = Path(self.config.get('project_path', '/tmp/default_project'))
        if not self.project_path.exists():
            self.project_path.mkdir(parents=True, exist_ok=True)

    def log_start(self):
        """Log the start of the workflow."""
        print(f"Starting workflow: {self.workflow_name}")

    def log_end(self):
        """Log the end of the workflow."""
        print(f"Completed workflow: {self.workflow_name}")

    def run(self):
        """Run the workflow. Must be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement the run method.")