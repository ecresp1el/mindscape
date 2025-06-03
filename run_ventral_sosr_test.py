import subprocess
import sys
import shutil
from pathlib import Path
from datetime import datetime as dt
import logging
logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger("ventral_sosr_test")

from mindscape.tools import add_workflow_to_config 




# New: Define test directory name to isolate all outputs of this script
test_turbo_subdir = "test_runs"  # Customize as needed
base_turbo_path = Path("/nfs/turbo/umms-parent") / test_turbo_subdir

### section: Setup paths and check for existing files ###
repo_root = Path(__file__).resolve().parent # Assuming this script is in the root of the repository
generated_runner_file = repo_root / f"run_workflows_VentralSOSRSTest.py" # Generated runner script path
generated_pipeline_file = repo_root / "mindscape" / "bioinformatics_workflow_engine" / "pipelines" / "testing_ventral_workflow.py"

# Full project output path
today = dt.today().strftime("%Y-%m-%d")
project_path = base_turbo_path / f"ventral_sosr_test-MannyCrespo-{today}"


# Track deleted paths
deleted_paths = []

# Check for existing test folder or generated files
existing_items = []
if project_path.exists():
    existing_items.append(f"  - {project_path}")
if generated_runner_file.exists():
    existing_items.append(f"  - {generated_runner_file}")
if generated_pipeline_file.exists():
    existing_items.append(f"  - {generated_pipeline_file}")

if existing_items:
    warning_msg = "⚠️  Detected existing test outputs in:\n" + "\n".join(existing_items)
    warning_msg += "\n❓ Delete these and continue? (y/N): "
    user_input = input(warning_msg).strip().lower()
    if user_input == "y":
        if project_path.exists():
            shutil.rmtree(project_path)
            logger.info(f"🗑️ Deleted test project directory: {project_path}")
            deleted_paths.append(str(project_path))
        if generated_runner_file.exists():
            generated_runner_file.unlink()
            logger.info(f"🗑️ Deleted runner script: {generated_runner_file}")
            deleted_paths.append(str(generated_runner_file))
        if generated_pipeline_file.exists():
            generated_pipeline_file.unlink()
            logger.info(f"🗑️ Deleted scaffolded pipeline file: {generated_pipeline_file}")
            deleted_paths.append(str(generated_pipeline_file))
        if deleted_paths:
            logger.info("🧾 Confirmed deleted paths:")
            for path in deleted_paths:
                logger.info(f"  - {path}")
    else:
        logger.info("❌ Aborting test to preserve existing files.")
        sys.exit(1)

# Script to automate project creation and runner execution using testscript_cli.py

def run_ventral_sosr_test():
    #cli_script = "testscript_cli.py"
    project_name = "ventral_sosr_test"
    experimenter_name = "MannyCrespo"
    email = "elcrespo@umich.edu"
    runner_name = "VentralSOSRSTest"

    command = [
        sys.executable, "testscript_cli.py",
        "--project_name", project_name,
        "--experimenter_name", experimenter_name,
        "--email", email,
        "--blank",
        "--blank_runner", runner_name,
        "--test_turbo_subdir", test_turbo_subdir
    ]

    logger.info("🚀 Launching automated project creation and workflow runner...")
    logger.debug("🔧 Command: " + " ".join(command))

    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        logger.info("✅ Project creation script executed successfully.")
        logger.debug("STDOUT:\n" + result.stdout)
        logger.debug("STDERR:\n" + result.stderr)
    except subprocess.CalledProcessError as e:
        logger.info(f"❌ Failed to execute testscript_cli.py: {e}")
        logger.debug("STDOUT:\n" + (e.stdout or ""))
        logger.debug("STDERR:\n" + (e.stderr or ""))
        return

    # Add existing CellRangerWorkflow to config
    try:
        logger.info("➕ Adding existing CellRangerWorkflow to config...")
        
        add_workflow_to_config.main([
            "--project_path", str(project_path),
            "--workflow_name", "CellRangerWorkflow"
        ])
        
        logger.info("✅ CellRangerWorkflow added successfully.")
    except Exception as e:
        logger.info(f"❌ Failed to add CellRangerWorkflow: {e}")
        return

    # Add new TestingVentralWorkflow to config
    try:
        logger.info("➕ Adding TestingVentralWorkflow to config...")
        
        add_workflow_to_config.main([
            "--project_path", str(project_path),
            "--workflow_name", "TestingVentralWorkflow"
        ])
        logger.info("✅ TestingVentralWorkflow added successfully.")
    except Exception as e:
        logger.info(f"❌ Failed to add TestingVentralWorkflow: {e}")
        return
    
    # Scaffold TestingVentralWorkflow (this generates a .py file dynamically and must match snake_case naming)
    scaffold_script = "mindscape/tools/create_workflow_scaffold.py"
    try:
        scaffold_cmd = [
            sys.executable, scaffold_script,
            "--name", "TestingVentralWorkflow"
        ]
        logger.info("🛠️  Scaffolding TestingVentralWorkflow...")
        result = subprocess.run(scaffold_cmd, check=True, text=True, capture_output=True)
        logger.info("✅ TestingVentralWorkflow scaffold created successfully.")
        logger.debug("STDOUT:\n" + result.stdout)
        logger.debug("STDERR:\n" + result.stderr)
    except subprocess.CalledProcessError as e:
        logger.info(f"❌ Failed to scaffold TestingVentralWorkflow: {e}")
        logger.debug("STDOUT:\n" + (e.stdout or ""))
        logger.debug("STDERR:\n" + (e.stderr or ""))
        return

    # Run the generated runner script (now includes --project_path as required by updated workflow base class)
    runner_script_path = Path("mindscape/bioinformatics_workflow_engine") / f"run_workflows_{runner_name}.py"
    run_runner_cmd = [
        sys.executable,
        str(runner_script_path),
        "--project_path",
        str(project_path)
    ]
    try:
        logger.info(f"🏃 Running the generated runner script: {runner_script_path} ...")
        result = subprocess.run(run_runner_cmd, check=True, text=True, capture_output=True)
        logger.info("✅ Runner script executed successfully.")
        logger.debug("STDOUT:\n" + result.stdout)
        logger.debug("STDERR:\n" + result.stderr)
    except subprocess.CalledProcessError as e:
        logger.info(f"❌ Failed to execute runner script {runner_script_path}: {e}")
        logger.debug("STDOUT:\n" + (e.stdout or ""))
        logger.debug("STDERR:\n" + (e.stderr or ""))
        return

if __name__ == "__main__":
    logger.info(f"🧹 Cleaning up previous test state...")
    run_ventral_sosr_test()
