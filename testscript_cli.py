import subprocess
import os
from datetime import datetime as dt
from pathlib import Path
from mindscape.config import create_config_template, write_config

def test_create_new_project_cli():
    """
    Test the MindScape CLI for creating a new project on Turbo storage.
    """
    # Define test inputs
    project_name = "TestProjectCLI"
    experimenter_name = "TestUser"
    working_directory = "/nfs/turbo/umms-parent/"  # Turbo directory

    # Ensure the working directory exists
    os.makedirs(working_directory, exist_ok=True)

    # Construct the CLI command
    command = [
        "python", "-m", "mindscape", "create-project",
        project_name, experimenter_name,
        "--working-directory", working_directory
    ]

    print(f"Running command: {' '.join(command)}")

    # Run the command and capture the output
    result = subprocess.run(command, capture_output=True, text=True)

    # Print the output for debugging
    print("STDOUT:")
    print(result.stdout)
    print("STDERR:")
    print(result.stderr)

    # Check if the command was successful
    if result.returncode == 0:
        print(f"✅ CLI test passed: Project '{project_name}' created successfully!")
    else:
        print(f"❌ CLI test failed with return code {result.returncode}.")
        print("Check the error messages above for details.")

    # Verify the project directory was created
    project_path = os.path.join(working_directory, f"{project_name}-{experimenter_name}-2025-05-29")
    if os.path.exists(project_path):
        print(f"✅ Project directory exists: {project_path}")
    else:
        print(f"❌ Project directory does not exist: {project_path}")

    # Debug: List the contents of the working directory
    print("\nContents of the working directory:")
    print(os.listdir(working_directory))


def create_new_project(
    project: str,
    experimenter: str,
    working_directory: str | None = None,
):
    """
    Creates a new MindScape project directory with the necessary structure.
    """
    # Debug: Print the working directory
    print(f"DEBUG: Received working_directory = {working_directory}")

    # Get the current date
    date = dt.today().strftime("%Y-%m-%d")

    # Set the working directory
    if working_directory is None:
        working_directory = "."
    wd = Path(working_directory).resolve()

    # Debug: Print the resolved working directory
    print(f"DEBUG: Resolved working_directory = {wd}")

    # Create the project name and path
    project_name = f"{project}-{experimenter}-{date}"
    project_path = wd / project_name

    # Debug: Print the project path
    print(f"DEBUG: Project path = {project_path}")

    # Check if the project already exists
    if project_path.exists():
        print(f"⚠️ Project '{project_name}' already exists at {project_path}.")
        return str(project_path)

    # Try creating the project directory
    try:
        print(f"DEBUG: Attempting to create project directory at {project_path}")
        project_path.mkdir(parents=True, exist_ok=True)
        (project_path / "data").mkdir(exist_ok=True)
        (project_path / "results").mkdir(exist_ok=True)
        (project_path / "logs").mkdir(exist_ok=True)
        (project_path / "config").mkdir(exist_ok=True)
        print(f"DEBUG: Successfully created project directory at {project_path}")
    except Exception as e:
        print(f"❌ Failed to create project directories: {e}")
        return None

    # Generate the configuration file
    config_path = project_path / "config/config.yaml"
    cfg, ruamelFile = create_config_template()
    cfg["project_name"] = project
    cfg["experimenter"] = experimenter
    cfg["date"] = date
    cfg["project_path"] = str(project_path)
    write_config(config_path, cfg)

    print(f"✅ New MindScape project created at: {project_path}")
    return str(project_path)


if __name__ == "__main__":
    test_create_new_project_cli()