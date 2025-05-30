import unittest
from utils.logger import setup_logger
from utils.paths import get_project_root

class TestUtils(unittest.TestCase):

    def setUp(self):
        # Provide a valid log file path for testing
        self.logger = setup_logger("test_logger", "test_log.log")
        self.project_root = get_project_root()

    def test_logger_setup(self):
        self.assertIsNotNone(self.logger)
        self.logger.info("Logger setup test passed.")

    def test_get_project_root(self):
        self.assertTrue(self.project_root.exists())
        self.logger.info(f"Project root is: {self.project_root}")

if __name__ == "__main__":
    unittest.main()