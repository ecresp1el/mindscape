
from mindscape.bioinformatics_workflow_engine.pipelines.my_test_workflow import FailingWorkflow
import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project_path", type=str, required=True)
    args = parser.parse_args()
    wf = FailingWorkflow(args.project_path)
    try:
        wf.run()
    except Exception as e:
        wf.mark_failed()
        print(f"[DEBUG] Exception in FailingWorkflow: {e}", file=sys.stderr)
        raise

if __name__ == "__main__":
    main()
