import unittest
import sys
from pathlib import Path

# Add the project root directory to sys.path
sys.path.insert(0, str(Path(__file__).parent.parent))

from pipelines.base_workflow import BaseWorkflow
import yaml

class TestBaseWorkflow(unittest.TestCase):

    def setUp(self):
        # Provide an absolute path to the config file
        self.config_path = Path(__file__).parent.parent / 'config/default_config.yaml'
        self.workflow = BaseWorkflow(config_path=self.config_path)

        # Load the configuration file to get the expected project_path
        with open(self.config_path, 'r') as file:
            self.config = yaml.safe_load(file)

    def test_load_config(self):
        # Test that the configuration is loaded as a dictionary
        self.assertIsInstance(self.workflow.config, dict)

    def test_set_up_paths(self):
        # Dynamically use the project_path from the configuration file
        expected_path = Path(self.config['project_path'])
        self.workflow.setup_paths()
        self.assertEqual(self.workflow.project_path, expected_path)
        self.assertTrue(self.workflow.project_path.exists())

    def test_run_method(self):
        # Test that the run method raises NotImplementedError
        with self.assertRaises(NotImplementedError):
            self.workflow.run()

if __name__ == '__main__':
    unittest.main()