import subprocess
import tempfile
from pathlib import Path
import shutil
import os

REPO_ROOT = Path(__file__).resolve().parents[3]

def test_baseworkflow_end_to_end():
    # Create a temporary test directory
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        project_name = "TestWorkflowProject"
        experimenter = "TestUser"
        runner_name = "MyTestWorkflow"

        testscript_cli_path = REPO_ROOT / "testscript_cli.py"
        print(f"Running project creation with script: {testscript_cli_path}")
        # Run project creation
        subprocess.run([
            "python", str(testscript_cli_path),
            "--project_name", project_name,
            "--experimenter_name", experimenter,
            "--blank",
            "--blank_runner", runner_name,
            "--mindscape_dry_run"
        ], check=True)

        # Locate created project path
        project_dir = next((p for p in Path("/nfs/turbo/umms-parent").glob(f"{project_name}-{experimenter}-*")), None)
        assert project_dir is not None, "Project directory was not created."
        print(f"Located project directory at: {project_dir}")

        add_workflow_script = REPO_ROOT / "mindscape" / "tools" / "add_workflow_to_config.py"
        print(f"Adding workflow to config using script: {add_workflow_script}")
        # Add workflow to config
        subprocess.run([
            "python", str(add_workflow_script),
            "--project_path", str(project_dir),
            "--workflow_name", runner_name
        ], check=True)

        runner_path = REPO_ROOT / "mindscape" / "bioinformatics_workflow_engine" / f"run_workflows_{runner_name}.py"
        print(f"Running generated runner script: {runner_path}")
        # Run the generated runner
        subprocess.run([
            "python", str(runner_path),
            "--project_path", str(project_dir)
        ], check=True)

        # Validate output
        logs_dir = project_dir / "logs"
        assert (logs_dir / f"{runner_name}.in_progress").exists(), ".in_progress file missing"
        assert (logs_dir / f"{runner_name}.completed").exists(), ".completed file missing"
        assert (logs_dir / "workflow_manager.log").exists(), "workflow_manager.log missing"

        log_text = (logs_dir / "workflow_manager.log").read_text()
        assert f"Starting workflow: {runner_name}" in log_text
        assert f"Completed workflow: {runner_name}" in log_text

        print("✅ BaseWorkflow integration test passed.")

if __name__ == "__main__":
    # Fallback: define a trivial MyTestWorkflow class if not already present
    try:
        from mindscape.bioinformatics_workflow_engine.baseworkflow import BaseWorkflow
    except ImportError:
        BaseWorkflow = object

    class MyTestWorkflow(BaseWorkflow):
        def run(self):
            print("Running fallback MyTestWorkflow.")
            return True

    print("🧪 Starting BaseWorkflow integration test...")
    try:
        test_baseworkflow_end_to_end()
        print("✅ test_baseworkflow_end_to_end() finished successfully.")
    except AssertionError as e:
        print("❌ Assertion failed during test_baseworkflow_end_to_end():", e)
        import traceback
        traceback.print_exc()
    except Exception as ex:
        print("❌ Exception occurred during test_baseworkflow_end_to_end():", ex)
        import traceback
        traceback.print_exc()