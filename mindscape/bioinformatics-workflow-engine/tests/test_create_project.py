import unittest
from pathlib import Path
from datetime import datetime as dt
from mindscape.create_project.new import create_new_project

class TestCreateProject(unittest.TestCase):

    def setUp(self):
        # Set up test parameters
        self.project_name = "TestProject"
        self.experimenter_name = "TestUser"
        self.turbo_shared_directory = "/tmp"  # Use /tmp for testing to avoid modifying real directories
        self.date = dt.today().strftime("%Y-%m-%d")  # Get the full date (YYYY-MM-DD)

    def test_create_new_project(self):
        # Call the create_new_project function
        project_path = create_new_project(
            project=self.project_name,
            experimenter=self.experimenter_name,
            working_directory=self.turbo_shared_directory,
        )

        # Verify the project path
        expected_path = Path(self.turbo_shared_directory) / f"{self.project_name}-{self.experimenter_name}-{self.date}"
        self.assertEqual(Path(project_path), expected_path)

        # Verify the directory structure
        self.assertTrue((Path(project_path) / "data").exists())
        self.assertTrue((Path(project_path) / "results").exists())
        self.assertTrue((Path(project_path) / "logs").exists())
        self.assertTrue((Path(project_path) / "config").exists())

    def tearDown(self):
        # Clean up the created project directory
        project_path = Path(self.turbo_shared_directory) / f"{self.project_name}-{self.experimenter_name}-{self.date}"
        if project_path.exists():
            for sub in project_path.iterdir():
                if sub.is_dir():
                    for file in sub.iterdir():
                        file.unlink()
                    sub.rmdir()
                else:
                    sub.unlink()
            project_path.rmdir()

if __name__ == "__main__":
    unittest.main()