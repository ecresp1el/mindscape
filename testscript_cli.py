"""This script demonstrates how to create a new MindScape project and run bioinformatics workflows.

Usage:
    python testscript_cli.py [--project_name PROJECT_NAME] [--experimenter_name EXPERIMENTER_NAME]

Description:
    - Creates a new project directory structure in a specified shared directory.
    - Initializes a configuration file for the project.
    - Runs the MindScape bioinformatics workflows on the created project.

Customization:
    - Modify 'project_name' and 'experimenter_name' to reflect your project details.
      You can also specify these as command line arguments:
        --project_name PROJECT_NAME
        --experimenter_name EXPERIMENTER_NAME
    - Set 'turbo_shared_directory' to the path where you want the project directory to be created.
    - Ensure that the MindScape module is installed and accessible in your Python environment.

Alternate CLI usage:
    You can also run the workflows with additional options using the command line interface:
    
    --email EMAIL_ADDRESS
        Specify an email address to receive notifications about the workflow status.
        This can be passed as a command line argument or set as an environment variable:
            export MINDSCAPE_EMAIL=your_email@example.com

    --mindscape_dry_run
        Perform a dry run of the workflows without executing them.
        This flag can be passed directly via the CLI or set as an environment variable:
            export MINDSCAPE_DRY_RUN=1

## 🧪 Example CLI Usage

Run the script with default project and experimenter names:
    python testscript_cli.py

Specify a custom project name and experimenter name:
    python testscript_cli.py --project_name MyProject --experimenter_name Alice

Include an email address for notifications:
    python testscript_cli.py --email alice@umich.edu

Perform a dry run without executing workflows:
    python testscript_cli.py --mindscape_dry_run

Combine all options:
    python testscript_cli.py --project_name MyProject --experimenter_name Alice --email alice@umich.edu --mindscape_dry_run

Set environment variables instead of passing arguments:
    export MINDSCAPE_EMAIL=alice@umich.edu
    export MINDSCAPE_DRY_RUN=1
    python testscript_cli.py --project_name MyProject --experimenter_name Alice
"""

import os, subprocess, sys
from datetime import datetime as dt
from pathlib import Path
import argparse
import json
import yaml

import mindscape as ms # Importing the main MindScape module

print("Imported MindScape!")

parser = argparse.ArgumentParser(description="Create a MindScape project and run workflows.")
parser.add_argument("--email", type=str, default=os.getenv("MINDSCAPE_EMAIL"), help="Email address for notifications")
parser.add_argument("--mindscape_dry_run", action="store_true", default=bool(int(os.getenv("MINDSCAPE_DRY_RUN", "0"))), help="Perform a dry run without executing workflows")
parser.add_argument("--project_name", type=str, default="TestProject", help="Name of the project")
parser.add_argument("--experimenter_name", type=str, default="TestUser", help="Name of the experimenter")
args = parser.parse_args()

print("Creating a new Mindscape project...")

# Define project and experimenter names from CLI arguments
project_name = args.project_name
experimenter_name = args.experimenter_name

# Define the directory where the project will be created (shared NFS directory)
turbo_shared_directory = "/nfs/turbo/umms-parent/"

## this will create a new project directory structure like:
# /nfs/turbo/umms-parent/TestProject-TestUser-2023-10-01/
# and initialize a configuration file in the specified working directory.

# Create a new MindScape project and get the path to the configuration file
path_config_file = ms.create_new_project(
    project=project_name,
    experimenter=experimenter_name,
    working_directory=turbo_shared_directory,
    email=args.email
)

print(f"Project created successfully at: {path_config_file}")
# Note: The above function creates a new project directory structure
# and initializes a configuration file in the specified working directory.
# The project is created in the specified NFS directory, which is shared
# across multiple nodes in the cluster.

print("Project creation completed successfully!")

# If mindscape_dry_run is set, update the config file's dry_run field before running the workflow
if args.mindscape_dry_run:
    config_path = Path(path_config_file) / "config/config.yaml"
    if not os.path.isfile(config_path):
        print(f"Error: Configuration file '{config_path}' does not exist. Cannot perform dry run update.")
    else:
        try:
            with open(config_path, "r") as f:
                config_data = yaml.safe_load(f)
        except Exception as read_err:
            print(f"Error reading configuration file '{config_path}': {read_err}")
        else:
            config_data["dry_run"] = True
            try:
                with open(config_path, "w") as f:
                    yaml.dump(config_data, f)
                print("Dry run enabled: updated configuration file successfully.")
            except Exception as write_err:
                print(f"Error writing to configuration file '{config_path}': {write_err}")

# Run the workflow script using subprocess to process the created project
try:
    subprocess.run(
        [
            "python",
            "mindscape/bioinformatics_workflow_engine/run_workflows.py",
            "--project_path",
            str(path_config_file)  # Pass the dynamically determined project path
        ],
        check=True
    )
    print("Workflows executed successfully!")
except subprocess.CalledProcessError as e:
    print(f"Error while running workflows: {e}")