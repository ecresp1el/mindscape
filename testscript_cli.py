"""This script demonstrates how to create a new MindScape project and run bioinformatics workflows.

Usage:
    python testscript_cli.py [OPTIONS]

Options:
    --project_name PROJECT_NAME
        Name of the project. Default: "TestProject"

    --experimenter_name EXPERIMENTER_NAME
        Name of the experimenter. Default: "TestUser"

    --email EMAIL_ADDRESS
        Email for SLURM notifications. Default: Value from $MINDSCAPE_EMAIL or None

    --mindscape_dry_run
        Perform a dry run (simulate workflow without execution).
        Default: False (can also set $MINDSCAPE_DRY_RUN=1)

    --blank
        Use a blank configuration (no default workflows included).
        Default: False

    --blank_runner [NAME]
        Generate a blank run_workflows_<NAME>.py runner template inside bioinformatics_workflow_engine/.
        If NAME is omitted, defaults to run_workflows_<project_name>_template.py.
        Useful for custom workflow execution scaffolds.

Description:
    ┌─➤ CLI Call: testscript_cli.py [with CLI args]
    │   ├── --project_name           # required
    │   ├── --experimenter_name      # required
    │   ├── --email                  # optional or via $MINDSCAPE_EMAIL
    │   ├── --mindscape_dry_run      # optional or via $MINDSCAPE_DRY_RUN
    │   ├── --blank                  # optional → use blank config
    │   └── --blank_runner           # optional → generate and run custom runner
    │
    ├── SELECT CONFIG FUNCTION:
    │   ├── If `--blank` → use create_blank_config_template
    │   └── else         → use create_config_template
    │
    ├── CALL: create_new_project(...)
    │   ├── Creates folder + config/config.yaml
    │   └── Returns: full project path
    │
    ├── IF --mindscape_dry_run or env var set:
    │   └── Rewrites config["dry_run"] = true and saves
    │
    └── WORKFLOW EXECUTION BRANCHES:
        |
        ├── IF `--blank_runner` is given:
        │   ├── Generate custom runner file:
        │   │   └── run_workflows_<custom_name>.py
        │   └── subprocess.run([
        │         "python",
        │         "mindscape/bioinformatics_workflow_engine/run_workflows_<custom_name>.py",
        │         "--project_path", <created_project_path>
        │       ])
        │
        └── ELSE (default case):
            └── subprocess.run([
                "python",
                "mindscape/bioinformatics_workflow_engine/run_workflows.py",
                "--project_path", <created_project_path>
                ])

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

Generate a blank project and custom runner:
    python testscript_cli.py --project_name MyBlankProject --experimenter_name Alice --blank --blank_runner MyRunner
"""

import os, subprocess, sys
from datetime import datetime as dt
from pathlib import Path
import argparse
import json
import yaml

import mindscape as ms # Importing the main MindScape module
from mindscape.utils.auxilliaryfunctions import create_blank_config_template, create_config_template
from mindscape.tools.generate_workflow_runner import generate_runner_template

print("Imported MindScape!")

def main():
    print("DEBUG: sys.argv =", sys.argv)
    parser = argparse.ArgumentParser(description="Create a MindScape project and run workflows.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--email", type=str, default=os.getenv("MINDSCAPE_EMAIL"), help="Email address for notifications")
    parser.add_argument("--mindscape_dry_run", action="store_true", default=bool(int(os.getenv("MINDSCAPE_DRY_RUN", "0"))), help="Perform a dry run without executing workflows")
    parser.add_argument("--project_name", type=str, default="TestProject", help="Name of the project")
    parser.add_argument("--experimenter_name", type=str, default="TestUser", help="Name of the experimenter")
    parser.add_argument("--blank", action="store_true", help="Create project with a blank workflow list")
    parser.add_argument("--blank_runner", nargs="?", const=True, help="Generate a run_workflows_template.py file. Optionally provide a custom name, e.g. --blank_runner MyRunner")
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


    # Select the appropriate config template function based on --blank flag
    config_function = create_blank_config_template if args.blank else create_config_template

    # Create a new MindScape project and get the path to the configuration file
    path_config_file = ms.create_new_project(
        project=project_name,
        experimenter=experimenter_name,
        working_directory=turbo_shared_directory,
        email=args.email,
        custom_config_function=config_function
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

    # Determine which workflow runner script to use
    if args.blank_runner:

        print("DEBUG: args.blank_runner =", args.blank_runner)

        runner_name = (
            f"run_workflows_{project_name.lower()}_template.py"
            if isinstance(args.blank_runner, bool)
            else f"run_workflows_{args.blank_runner}.py"
        )
        print("DEBUG: Generated runner_name =", runner_name)

        runner_script_path = Path("mindscape/bioinformatics_workflow_engine") / runner_name
        print("DEBUG: runner_script_path =", runner_script_path)

        generate_runner_template(output_path=runner_script_path)

        runner_script_absolute = Path.cwd() / runner_script_path
        print(f"✅ Template '{runner_script_path.name}' created. You can now add your workflows to it.")
        print("📂 Executing runner at:", runner_script_absolute)

        try:
            with open(runner_script_absolute, "r") as f:
                preview_lines = "".join(f.readlines()[:5])
                print("🧾 Preview runner script:\n", preview_lines)
        except Exception as e:
            print(f"⚠️ Failed to preview runner: {e}")

        print("DEBUG: Launching subprocess for", runner_script_absolute)
        result = subprocess.run(
            [
                sys.executable,
                str(runner_script_absolute),
                "--project_path",
                str(path_config_file)
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        print("DEBUG: Runner STDOUT:")
        print(result.stdout)
        print("DEBUG: Runner STDERR:")
        print(result.stderr)
        if result.returncode != 0:
            print(f"❌ Error: Runner exited with return code {result.returncode}")
            sys.exit(result.returncode)
        else:
            sys.exit(0)  # ✅ Prevent fallthrough into default runner
    else:
        runner_script_path = Path("mindscape/bioinformatics_workflow_engine/run_workflows.py")

        print("DEBUG: Using default runner =", runner_script_path)

        try:
            subprocess.run(
                [
                    sys.executable,
                    os.path.abspath(str(runner_script_path)),
                    "--project_path",
                    str(path_config_file)
                ],
                check=True
            )
            print("Workflows executed successfully!")
        except subprocess.CalledProcessError as e:
            print(f"Error while running workflows: {e}")

if __name__ == "__main__":
    main()