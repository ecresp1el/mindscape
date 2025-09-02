


import warnings

def warn_if_missing_from_config(pipeline_dir, configured_names):
    """
    Emits a warning if a Python workflow module is not present in the config.

    Args:
        pipeline_dir (Path): Directory containing all pipeline workflow scripts.
        configured_names (list[str]): List of workflow names configured in the config file.

    Example:
        If `pipeline_dir` contains `cell_ranger_workflow.py` and `qc_workflow.py`
        but the config only includes `CellRangerWorkflow`, this function will print
        a warning for `QCWorkflow`.

    This helps users detect workflows they might have forgotten to include.
    """
    from pathlib import Path

    # Find all Python files ending in "_workflow.py"
    all_workflows = {
        f.stem.replace("_workflow", "").title().replace("_", "") + "Workflow"
        for f in Path(pipeline_dir).glob("*_workflow.py")
    }

    # Remove duplicates from config and normalize case
    config_set = {str(name).strip() for name in configured_names}
    missing = all_workflows - config_set

    if missing:
        warnings.warn(
            f"The following workflow(s) exist in the pipeline directory but are not enabled in the config:\n  - "
            + "\n  - ".join(sorted(missing))
        )