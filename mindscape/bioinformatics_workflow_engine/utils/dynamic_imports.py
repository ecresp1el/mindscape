import importlib
import inspect
from pathlib import Path
import sys
from mindscape.bioinformatics_workflow_engine.pipelines.base_workflow import BaseWorkflow

def dynamic_import_workflows(pipelines_dir: Path):
    workflows = {}
    for py_file in pipelines_dir.glob("*.py"):
        if py_file.name.startswith("__") or py_file.name == "base_workflow.py":
            continue

        module_path = f"pipelines.{py_file.stem}"
        try:
            module = importlib.import_module(module_path)
            for name, obj in inspect.getmembers(module, inspect.isclass):
                if issubclass(obj, BaseWorkflow) and obj.__name__ != "BaseWorkflow":
                    workflows[name] = obj
        except Exception as e:
            print(f"⚠️ Could not import {module_path}: {e}")
    return workflows