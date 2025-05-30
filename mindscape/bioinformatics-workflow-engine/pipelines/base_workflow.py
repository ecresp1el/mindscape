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
        if not self.config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.config_path}")
        with open(self.config_path, 'r') as file:
            return yaml.safe_load(file)

    def setup_logging(self):
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(self.__class__.__name__)

    def run(self):
        raise NotImplementedError("Subclasses must implement the run method.")