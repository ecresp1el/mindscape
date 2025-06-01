from __future__ import annotations

import os
import warnings
from pathlib import Path
from ruamel.yaml import YAML
import yaml

def create_config_template(project_path: str, project_name: str, experimenter: str, date: str, email: str ="elcrespo@umich.edu"):
    """
    Creates a complete and usable config.yaml with default values (SLURM disabled).
    """
    yaml_str = f"""\
project_name: {project_name}
experimenter: {experimenter}
date: {date}

project_path: {project_path}

data_dir: data
results_dir: results
logs_dir: logs


use_slurm: false
dry_run: false
force_rerun: false
email: {email}

slurm:
  cpus: "8"
  mem: "64G"
  time: "12:00:00"

parameters:
  analysis_type: default
  threshold: 0.5
  max_iterations: 1000

workflows:
  - name: CellRangerWorkflow
    enabled: true
  - name: QCWorkflow
    enabled: true
  - name: VentralWorkflow
    enabled: false
"""
    ruamelFile = YAML()
    ruamelFile.preserve_quotes = True
    cfg_file = ruamelFile.load(yaml_str)
    return cfg_file, ruamelFile


def read_config(configname: str | Path) -> dict:
    """
    Reads a structured config file defining a project.

    Parameters
    ----------
    configname : str or Path
        Path to the configuration file.

    Returns
    -------
    dict
        Parsed configuration as a dictionary.
    """
    ruamelFile = YAML()
    path = Path(configname)
    if not path.exists():
        raise FileNotFoundError(
            f"Config file at {path} not found. Please make sure the file exists and/or the path is correct!"
        )

    with open(path, "r") as f:
        cfg = ruamelFile.load(f)

    # Ensure the project path matches the current directory
    curr_dir = str(path.parent.resolve())
    if cfg.get("project_path") != curr_dir:
        cfg["project_path"] = curr_dir
        write_config(configname, cfg)

    return cfg


def write_config(configname: str | Path, cfg: dict):
    """
    Writes a structured config file to disk.

    Parameters
    ----------
    configname : str or Path
        Path to the configuration file.
    cfg : dict
        Fully formed configuration data to write.
    """
    yaml = YAML()
    yaml.preserve_quotes = True
    with open(configname, "w") as f:
        yaml.dump(cfg, f)


def edit_config(configname: str | Path, edits: dict, output_name: str = "") -> dict:
    """
    Edits and saves a config file from a dictionary.

    Parameters
    ----------
    configname : str or Path
        Path to the configuration file.
    edits : dict
        Key-value pairs to edit in the config.
    output_name : str, optional
        If provided, saves the edited config to a new file. Defaults to overwriting the original file.

    Returns
    -------
    dict
        Updated configuration.
    """
    cfg = read_config(configname)
    for key, value in edits.items():
        cfg[key] = value

    if not output_name:
        output_name = configname

    try:
        write_config(output_name, cfg)
    except Exception as e:
        warnings.warn(
            f"Some edits could not be written due to: {e}. The configuration file will remain unchanged."
        )
        for key in edits:
            cfg.pop(key)
        write_config(output_name, cfg)

    return cfg