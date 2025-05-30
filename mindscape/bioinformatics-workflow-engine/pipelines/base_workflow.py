from __future__ import annotations

import os
import logging
from pathlib import Path
import yaml

class BaseWorkflow:
    def __init__(self, config_path):
        self.config_path = config_path
        self.config = self.load_config()

    def load_config(self):
        """Load the YAML configuration file."""
        with open(self.config_path, 'r') as file:
            return yaml.safe_load(file)

    def setup_paths(self):
        """Set up project paths."""
        self.project_path = Path(self.config.get('project_path', '/tmp/project'))
        if not self.project_path.exists():
            self.project_path.mkdir(parents=True, exist_ok=True)

    def run(self):
        """Run the workflow. Must be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement the run method.")