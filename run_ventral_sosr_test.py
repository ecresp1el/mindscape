import subprocess
import sys
import shutil
from pathlib import Path
from datetime import datetime as dt

# New: Define test directory name to isolate all outputs
test_turbo_subdir = "test_runs"  # Customize as needed
base_turbo_path = Path("/nfs/turbo/umms-parent") / test_turbo_subdir
repo_root = Path(__file__).resolve().parent
generated_runner_file = repo_root / f"run_workflows_VentralSOSRSTest.py"
generated_pipeline_file = repo_root / "mindscape" / "bioinformatics_workflow_engine" / "pipelines" / "testing_ventral_workflow.py"

# Full project output path
today = dt.today().strftime("%Y-%m-%d")
project_path = base_turbo_path / f"ventral_sosr_test-MannyCrespo-{today}"

# Check for existing test folder or generated files
if project_path.exists() or generated_runner_file.exists() or generated_pipeline_file.exists():
    user_input = input(
        f"⚠️  Detected existing test outputs in:\n"
        f"  - {project_path if project_path.exists() else ''}\n"
        f"  - {generated_runner_file if generated_runner_file.exists() else ''}\n"
        f"  - {generated_pipeline_file if generated_pipeline_file.exists() else ''}\n"
        f"❓ Delete these and continue? (y/N): "
    ).strip().lower()
    if user_input == "y":
        if project_path.exists():
            shutil.rmtree(project_path)
            print(f"🗑️ Deleted test project directory: {project_path}")
        if generated_runner_file.exists():
            generated_runner_file.unlink()
            print(f"🗑️ Deleted runner script: {generated_runner_file}")
        if generated_pipeline_file.exists():
            generated_pipeline_file.unlink()
            print(f"🗑️ Deleted scaffolded pipeline file: {generated_pipeline_file}")
    else:
        print("❌ Aborting test to preserve existing files.")
        sys.exit(1)

# Script to automate project creation and runner execution using testscript_cli.py

def run_ventral_sosr_test():
    cli_script = "testscript_cli.py"
    project_name = "ventral_sosr_test"
    experimenter_name = "MannyCrespo"
    email = "elcrespo@umich.edu"
    runner_name = "VentralSOSRSTest"

    command = [
        sys.executable, cli_script,
        "--project_name", project_name,
        "--experimenter_name", experimenter_name,
        "--email", email,
        "--blank",
        "--blank_runner", runner_name,
        "--test_turbo_subdir", test_turbo_subdir
    ]

    print("🚀 Launching automated project creation and workflow runner...")
    print("🔧 Command:", " ".join(command))

    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        print("✅ Project creation script executed successfully.")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to execute testscript_cli.py: {e}")
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)
        return

    # Add existing CellRangerWorkflow to config
    add_workflow_script = "mindscape/tools/add_workflow_to_config.py"
    try:
        add_cellranger_cmd = [
            sys.executable, add_workflow_script,
            "--project_path", str(project_path),
            "--workflow_name", "CellRangerWorkflow"
        ]
        print("➕ Adding existing CellRangerWorkflow to config...")
        result = subprocess.run(add_cellranger_cmd, check=True, text=True, capture_output=True)
        print("✅ CellRangerWorkflow added successfully.")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to add CellRangerWorkflow: {e}")
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)
        return

    # Add new TestingVentralWorkflow to config
    try:
        add_ventral_cmd = [
            sys.executable, add_workflow_script,
            "--project_path", str(project_path),
            "--workflow_name", "TestingVentralWorkflow"
        ]
        print("➕ Adding new TestingVentralWorkflow to config...")
        result = subprocess.run(add_ventral_cmd, check=True, text=True, capture_output=True)
        print("✅ TestingVentralWorkflow added successfully.")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to add TestingVentralWorkflow: {e}")
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)
        return

    # Scaffold TestingVentralWorkflow
    scaffold_script = "mindscape/tools/create_workflow_scaffold.py"
    try:
        scaffold_cmd = [
            sys.executable, scaffold_script,
            "--name", "TestingVentralWorkflow"
        ]
        print("🛠️  Scaffolding TestingVentralWorkflow...")
        result = subprocess.run(scaffold_cmd, check=True, text=True, capture_output=True)
        print("✅ TestingVentralWorkflow scaffold created successfully.")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to scaffold TestingVentralWorkflow: {e}")
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)
        return

    # Run the generated runner script after workflows are added
    runner_script_path = Path("mindscape/bioinformatics_workflow_engine") / f"run_workflows_{runner_name}.py"
    run_runner_cmd = [sys.executable, str(runner_script_path)]
    try:
        print(f"🏃 Running the generated runner script: {runner_script_path} ...")
        result = subprocess.run(run_runner_cmd, check=True, text=True, capture_output=True)
        print("✅ Runner script executed successfully.")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to execute runner script {runner_script_path}: {e}")
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)
        return

if __name__ == "__main__":
    run_ventral_sosr_test()
