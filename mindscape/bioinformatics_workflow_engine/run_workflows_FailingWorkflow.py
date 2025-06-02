
from mindscape.bioinformatics_workflow_engine.pipelines.my_test_workflow import FailingWorkflow
import argparse
import sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_path", type=str, required=True)
    args = parser.parse_args()
    config_path = Path(args.project_path) / "config" / "config.yaml"
    wf = FailingWorkflow(config_path=config_path)
    try:
        wf.run()
    except Exception as e:
        # Use BaseWorkflow logic to get the marker path
        marker_path = wf.get_completion_marker_path()
        marker_path.parent.mkdir(parents=True, exist_ok=True)
        marker_path.touch()
        print(f"[DEBUG] Exception in FailingWorkflow: {e}", file=sys.stderr)
        print(f"[DEBUG] Wrote .failed marker to: {marker_path}")
        raise

if __name__ == "__main__":
    main()
