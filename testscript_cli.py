"""This script demonstrates how to create a new MindScape project and run bioinformatics workflows.

Usage:
    python testscript_cli.py

Description:
    - Creates a new project directory structure in a specified shared directory.
    - Initializes a configuration file for the project.
    - Runs the MindScape bioinformatics workflows on the created project.

Customization:
    - Modify 'project_name' and 'experimenter_name' to reflect your project details.
    - Set 'turbo_shared_directory' to the path where you want the project directory to be created.
    - Ensure that the MindScape module is installed and accessible in your Python environment.

Alternate CLI usage:
    You can also run the workflows with additional options using the command line interface:
    
    --email EMAIL_ADDRESS
        Specify an email address to receive notifications about the workflow status.
        This can be passed as a command line argument or set as an environment variable:
            export MIND_EMAIL=your_email@example.com

    --dry_run
        Perform a dry run of the workflows without executing them.
        This flag can be passed directly via the CLI or set as an environment variable:
            export MIND_DRY_RUN=1
"""

import os, subprocess, sys
from datetime import datetime as dt
from pathlib import Path
import argparse
import json

import mindscape as ms # Importing the main MindScape module

print("Imported MindScape!")

parser = argparse.ArgumentParser(description="Create a MindScape project and run workflows.")
parser.add_argument("--email", type=str, default=os.getenv("MIND_EMAIL"), help="Email address for notifications")
parser.add_argument("--dry_run", action="store_true", default=bool(int(os.getenv("MIND_DRY_RUN", "0"))), help="Perform a dry run without executing workflows")
args = parser.parse_args()

print("Creating a new project...")

# Define project and experimenter names
project_name = "TestProject"
experimenter_name = "TestUser"

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

# If dry_run is set, update the config file's dry_run field before running the workflow
if args.dry_run:
    try:
        with open(path_config_file, "r") as f:
            config_data = json.load(f)
        config_data["dry_run"] = True
        with open(path_config_file, "w") as f:
            json.dump(config_data, f, indent=4)
        print("Updated configuration file for dry run.")
    except Exception as e:
        print(f"Failed to update config file for dry run: {e}")

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