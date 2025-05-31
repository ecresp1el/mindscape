import os
import shutil
import warnings
from pathlib import Path
from datetime import datetime as dt
from mindscape.utils.auxilliaryfunctions import create_config_template, write_config

def create_new_project(
    project: str,
    experimenter: str,
    working_directory: str | None = None,
) -> str:
    """
    Creates a new MindScape project directory with the necessary structure.

    Args:
        project (str): The name of the project.
        experimenter (str): The name of the experimenter.
        working_directory (str | None): The base directory where the project will be created.

    Returns:
        str: The full path to the created project directory.
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
        #if it exists, we will not overwrite it
        warnings.warn(f"Project '{project_name}' already exists. No changes made.")
        return str(project_path)
    else:
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

