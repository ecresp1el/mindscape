import subprocess
import sys

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
        "--blank_runner", runner_name
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

    # Construct the full resolved project path
    from datetime import datetime as dt
    today = dt.today().strftime("%Y-%m-%d")
    project_path = f"/nfs/turbo/umms-parent/{project_name}-{experimenter_name}-{today}"

    # Add existing CellRangerWorkflow to config
    add_workflow_script = "mindscape/tools/add_workflow_to_config.py"
    try:
        add_cellranger_cmd = [
            sys.executable, add_workflow_script,
            "--project_path", project_path,
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

    # Add new VentralAnalysisWorkflow to config
    try:
        add_ventral_cmd = [
            sys.executable, add_workflow_script,
            "--project_path", project_path,
            "--workflow_name", "VentralAnalysisWorkflow"
        ]
        print("➕ Adding new VentralAnalysisWorkflow to config...")
        result = subprocess.run(add_ventral_cmd, check=True, text=True, capture_output=True)
        print("✅ VentralAnalysisWorkflow added successfully.")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to add VentralAnalysisWorkflow: {e}")
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)
        return

    # Scaffold VentralAnalysisWorkflow
    scaffold_script = "mindscape/tools/create_workflow_scaffold.py"
    try:
        scaffold_cmd = [
            sys.executable, scaffold_script,
            "--project_name", project_name,
            "--workflow_name", "VentralAnalysisWorkflow"
        ]
        print("🛠️  Scaffolding VentralAnalysisWorkflow...")
        result = subprocess.run(scaffold_cmd, check=True, text=True, capture_output=True)
        print("✅ VentralAnalysisWorkflow scaffold created successfully.")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to scaffold VentralAnalysisWorkflow: {e}")
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)
        return

    # Run the generated runner script after workflows are added
    runner_script = f"{runner_name}.py"
    run_runner_cmd = [sys.executable, runner_script]
    try:
        print(f"🏃 Running the generated runner script: {runner_script} ...")
        result = subprocess.run(run_runner_cmd, check=True, text=True, capture_output=True)
        print("✅ Runner script executed successfully.")
        print("STDOUT:\n", result.stdout)
        print("STDERR:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to execute runner script {runner_script}: {e}")
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)
        return

if __name__ == "__main__":
    run_ventral_sosr_test()
