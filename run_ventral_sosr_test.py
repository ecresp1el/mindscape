import subprocess
import sys
import shutil
from pathlib import Path
from datetime import datetime as dt
import logging
logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger("ventral_sosr_test")

from mindscape.tools import add_workflow_to_config, create_workflow_scaffold

# New: Define test directory name to isolate all outputs of this script
test_turbo_subdir = "test_runs"  # Customize as needed
base_turbo_path = Path("/nfs/turbo/umms-parent") / test_turbo_subdir

### section: Setup paths and check for existing files ###
repo_root = Path(__file__).resolve().parent # Assuming this script is in the root of the repository
generated_runner_file = repo_root / "mindscape" / "bioinformatics_workflow_engine" / f"run_workflows_VentralSOSRSTest.py"
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
    import os
    os.environ["PYTHONPATH"] = str(repo_root)

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
        logger.info("STDOUT:\n" + (e.stdout or "[no stdout]"))
        logger.info("STDERR:\n" + (e.stderr or "[no stderr]"))
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
    try:
        logger.info("🛠️  Scaffolding TestingVentralWorkflow...")
        create_workflow_scaffold.main([
            "--name", "TestingVentralWorkflow"
        ])
        logger.info("✅ TestingVentralWorkflow scaffold created successfully.")
    except Exception as e:
        logger.info(f"❌ Failed to scaffold TestingVentralWorkflow: {e}")
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

    # Preview loaded workflow classes and verify inheritance
    try:
        from mindscape.bioinformatics_workflow_engine.utils.dynamic_imports import dynamic_import_workflows
        from mindscape.bioinformatics_workflow_engine.pipelines.base_workflow import BaseWorkflow
        import yaml

        config_path = project_path / "config" / "config.yaml"
        with open(config_path, "r") as f:
            config_data = yaml.safe_load(f)

        workflow_names = config_data.get("workflows", [])
        pipelines_dir = repo_root / "mindscape" / "bioinformatics_workflow_engine" / "pipelines"
        
        available_workflows = dynamic_import_workflows(
            pipelines_dir,
            module_prefix="mindscape.bioinformatics_workflow_engine.pipelines"
        )
        logger.info("🔍 Previewing loaded workflow classes from config.yaml:")
        for wf in workflow_names:
            name = wf["name"] if isinstance(wf, dict) else wf
            wf_cls = available_workflows.get(name)
            if wf_cls is None:
                logger.warning(f"  ⚠️  Workflow '{name}' not found.")
            else:
                logger.info(f"  ✅ Loaded workflow class: {wf_cls.__name__} from {wf_cls.__module__}")
                if issubclass(wf_cls, BaseWorkflow):
                    logger.info(f"     🧬 Confirmed {wf_cls.__name__} inherits from BaseWorkflow")
                else:
                    logger.warning(f"     ❌ {wf_cls.__name__} does NOT inherit from BaseWorkflow")
    except Exception as e:
        logger.warning(f"⚠️ Could not verify workflow class inheritance: {e}")

    # Final cleanup prompt
    cleanup_prompt = input("\n🧪 This was a test run. Delete all generated files? (y/N): ").strip().lower()
    if cleanup_prompt == "y":
        paths_to_delete = []
        if project_path.exists():
            shutil.rmtree(project_path)
            paths_to_delete.append(project_path)
        if generated_runner_file.exists():
            generated_runner_file.unlink()
            paths_to_delete.append(generated_runner_file)
        if generated_pipeline_file.exists():
            generated_pipeline_file.unlink()
            paths_to_delete.append(generated_pipeline_file)

        if paths_to_delete:
            logger.info("🧼 Deleted the following test artifacts:")
            for path in paths_to_delete:
                if path == project_path:
                    logger.info(f"  - Project output directory: {path}")
                elif path == generated_runner_file:
                    logger.info(f"  - Generated runner script: {path}")
                elif path == generated_pipeline_file:
                    logger.info(f"  - Scaffolded workflow file: {path}")
                else:
                    logger.info(f"  - {path}")
    else:
        logger.info("🗂️ Retained test output files.")

if __name__ == "__main__":
    logger.info(f"🧹 Cleaning up previous test state...")
    run_ventral_sosr_test()
