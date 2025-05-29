import subprocess
import os

def test_create_new_project_cli():
    """
    Test the MindScape CLI for creating a new project.
    """
    # Define test inputs
    project_name = "TestProjectCLI"
    experimenter_name = "TestUser"
    working_directory = "/tmp/mindscape_test_cli"

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
    project_path = os.path.join(working_directory, f"{project_name}-{experimenter_name}-{os.getenv('DATE', '2025-05-29')}")
    if os.path.exists(project_path):
        print(f"✅ Project directory exists: {project_path}")
    else:
        print(f"❌ Project directory does not exist: {project_path}")


if __name__ == "__main__":
    test_create_new_project_cli()