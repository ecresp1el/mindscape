import unittest
from pipelines.base_workflow import BaseWorkflow
from pathlib import Path

class TestBaseWorkflow(unittest.TestCase):

    def setUp(self):
        # Provide an absolute path to the config file
        config_path = Path(__file__).parent.parent / 'config/default_config.yaml'
        self.workflow = BaseWorkflow(config_path=config_path)

    def test_load_config(self):
        # Test that the configuration is loaded as a dictionary
        self.assertIsInstance(self.workflow.config, dict)

    def test_set_up_paths(self):
        # Test setting up paths
        test_path = '/tmp/test_project_path'
        self.workflow.setup_paths()
        self.assertEqual(self.workflow.project_path, Path(test_path))
        self.assertTrue(self.workflow.project_path.exists())

    def test_run_method(self):
        # Test that the run method raises NotImplementedError
        with self.assertRaises(NotImplementedError):
            self.workflow.run()

if __name__ == '__main__':
    unittest.main()