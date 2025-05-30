from __future__ import annotations

import os
import logging
from pathlib import Path

class BaseWorkflow:
    def __init__(self, config_path):
        self.config_path = config_path
        self.config = self.load_config()

    def load_config(self):
        """Load the YAML configuration file."""
        import yaml
        with open(self.config_path, 'r') as file:
            return yaml.safe_load(file)

    def setup_paths(self):
        """Set up project paths."""
        # Use the project_path from the config or default to a writable directory
        self.project_path = Path(self.config.get('project_path', '/tmp/default_project'))
        if not self.project_path.exists():
            self.project_path.mkdir(parents=True, exist_ok=True)

    def run(self):
        """Run the workflow. Must be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement the run method.")