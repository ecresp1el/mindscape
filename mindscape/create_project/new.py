import os
import shutil
import warnings
from pathlib import Path
from datetime import datetime as dt
from mindscape.utils.auxilliaryfunction import create_config_template, write_config, read_config


def create_new_project(
    project: str,
    experimenter: str,
    working_directory: str | None = None,
):
    """
    Creates a new MindScape project directory with the necessary structure.

    Parameters
    ----------
    project : str
        Name of the project.
    experimenter : str
        Name of the experimenter.
    working_directory : str, optional
        Directory where the project will be created. Defaults to the current working directory.

    Returns
    -------
    str
        Path to the new project directory.
    """
    # Get the current date
    date = dt.today().strftime("%Y-%m-%d")

    # Set the working directory
    if working_directory is None:
        working_directory = "."
    wd = Path(working_directory).resolve()

    # Create the project name and path
    project_name = f"{project}-{experimenter}-{date}"
    project_path = wd / project_name

    # Check if the project already exists
    if project_path.exists():
        print(f"⚠️ Project '{project_name}' already exists at {project_path}.")
        return str(project_path)

    # Create project and subdirectories
    project_path.mkdir(parents=True, exist_ok=True)
    (project_path / "data").mkdir(exist_ok=True)
    (project_path / "results").mkdir(exist_ok=True)
    (project_path / "logs").mkdir(exist_ok=True)
    (project_path / "config").mkdir(exist_ok=True)

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

