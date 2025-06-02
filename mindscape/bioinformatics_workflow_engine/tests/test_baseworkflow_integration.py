import subprocess
import tempfile
from pathlib import Path
import shutil
import os

def test_baseworkflow_end_to_end():
    # Create a temporary test directory
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        project_name = "TestWorkflowProject"
        experimenter = "TestUser"
        runner_name = "MyTestWorkflow"

        # Run project creation
        subprocess.run([
            "python", "testscript_cli.py",
            "--project_name", project_name,
            "--experimenter_name", experimenter,
            "--blank",
            "--blank_runner", runner_name,
            "--mindscape_dry_run"
        ], check=True)

        # Locate created project path
        project_dir = next((p for p in Path("/nfs/turbo/umms-parent").glob(f"{project_name}-{experimenter}-*")), None)
        assert project_dir is not None, "Project directory was not created."

        # Add workflow to config
        subprocess.run([
            "python", "mindscape/tools/add_workflow_to_config.py",
            "--project_path", str(project_dir),
            "--workflow_name", runner_name
        ], check=True)

        # Run the generated runner
        runner_path = Path("mindscape/bioinformatics_workflow_engine") / f"run_workflows_{runner_name}.py"
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