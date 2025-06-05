# workflow/scripts/create_project.py
import os
import yaml
from pathlib import Path
from datetime import datetime

def main(config_path):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    date = datetime.today().strftime("%Y-%m-%d")
    project_name = config["project_name"]
    experimenter = config["experimenter"]
    turbo = Path(config["shared_turbo_directory"])
    full_name = f"{project_name}-{experimenter}-{date}"
    project_path = turbo / full_name

    # Create directory structure
    for subdir in ["data", "results", "logs", "config"]:
        (project_path / subdir).mkdir(parents=True, exist_ok=True)

    # Save resolved project_path into the same config for later rules
    config["project_path"] = str(project_path)
    with open(project_path / "config/config.yaml", "w") as out:
        yaml.dump(config, out)

    # ALSO: Write updated config back to the original input path
    with open(config_path, "w") as out_main:
        yaml.dump(config, out_main)

    print(f"✅ Created project structure at: {project_path}")
    
    # Create marker file for Snakemake completion
    
    #OLD: writes marker file to the github repo
    #marker_file = Path("results/create_project.done")
    # ✅ NEW: write marker to the *actual* project_path
    marker_file = project_path / "results" / "create_project.done"
    marker_file.parent.mkdir(parents=True, exist_ok=True)
    marker_file.touch()

# Snakemake provides a special object `snakemake` for scripts
main(snakemake.input[0])