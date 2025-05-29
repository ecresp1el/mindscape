import os, subprocess, sys
from datetime import datetime as dt
from pathlib import Path

import mindscape as ms # Importing the main MindScape module

print("Imported MindScape!")

print("Creating a new project...")

project_name = "TestProject"
experimenter_name = "TestUser"
turbo_shared_directory = "nfs/turbo/umms-parent/"

path_config_file = ms.create_new_project(
    project=project_name,
    experimenter=experimenter_name,
    working_directory=turbo_shared_directory
)

print(f"Project created successfully at: {path_config_file}")
# Note: The above function creates a new project directory structure
# and initializes a configuration file in the specified working directory.
# The project is created in the specified NFS directory, which is shared
# across multiple nodes in the cluster.
print("Project creation completed successfully!")