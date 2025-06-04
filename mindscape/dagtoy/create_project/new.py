# mindscape/dagtoy/create_project/new.py

from pathlib import Path
from datetime import datetime
import yaml


def create_config_template():
    return {
        "project_name": "TestProject",
        "experimenter": "Manny",
        "date": datetime.today().strftime("%Y-%m-%d"),
        "workflows": [
            {"name": "DataImportWorkflow", "enabled": True, "subworkflows": [], "depends_on": []},
            {"name": "AlignmentAndMoleculeCountingWorkflow", "enabled": True, "subworkflows": [], "depends_on": ["DataImportWorkflow"]},
            {"name": "CellFilteringWorkflow", "enabled": True, "subworkflows": [], "depends_on": ["AlignmentAndMoleculeCountingWorkflow"]},
            {"name": "DoubletScoringWorkflow", "enabled": True, "subworkflows": [], "depends_on": ["CellFilteringWorkflow"]},
            {"name": "CellSizeEstimationWorkflow", "enabled": True, "subworkflows": [], "depends_on": ["CellFilteringWorkflow"]},
            {"name": "GeneVarianceAnalysisWorkflow", "enabled": True, "subworkflows": [], "depends_on": ["CellFilteringWorkflow"]},
            {"name": "DimensionalityReductionWorkflow", "enabled": False, "subworkflows": [], "depends_on": ["GeneVarianceAnalysisWorkflow"]},
            {"name": "ManifoldRepresentationWorkflow", "enabled": True, "subworkflows": [], "depends_on": ["DimensionalityReductionWorkflow"]},
            {"name": "ClusteringAndDEWorkflow", "enabled": True, "subworkflows": [], "depends_on": ["ManifoldRepresentationWorkflow"]},
            {"name": "TrajectoryInferenceWorkflow", "enabled": True, "subworkflows": [], "depends_on": ["ClusteringAndDEWorkflow"]},
            {"name": "VelocityEstimationWorkflow", "enabled": True, "subworkflows": [], "depends_on": ["TrajectoryInferenceWorkflow"]},
            {"name": "CellTypeAnnotationWorkflow", "enabled": True, "subworkflows": [], "depends_on": ["ClusteringAndDEWorkflow"]},
            {"name": "IntegrationWorkflow", "enabled": True, "subworkflows": [], "depends_on": ["CellTypeAnnotationWorkflow"]},
            {"name": "MultiOmicsIntegrationWorkflow", "enabled": True, "subworkflows": [], "depends_on": ["IntegrationWorkflow"]}
        ],
        "slurm": {
            "cpus": 8,
            "mem": "32G",
            "time": "08:00:00"
        },
        "email": "elcrespo@umich.edu"
    }


def write_config(config_path: Path, config_dict: dict):
    config_path.parent.mkdir(parents=True, exist_ok=True)
    with open(config_path, "w") as f:
        yaml.dump(config_dict, f)


def create_new_project(base_dir="mindscape/dagtoy", project_name="TestProject", experimenter="Manny") -> tuple[Path, Path]:
    date_str = datetime.today().strftime("%Y-%m-%d")
    project_dir = Path(base_dir) / f"{project_name}-{experimenter}-{date_str}"

    # Create core directories
    for subdir in ["data", "results", "logs", "config"]:
        (project_dir / subdir).mkdir(parents=True, exist_ok=True)

    # Create config.yaml with defaults
    config_path = project_dir / "config" / "config.yaml"
    cfg = create_config_template()
    cfg["project_name"] = project_name
    cfg["experimenter"] = experimenter
    cfg["project_path"] = str(project_dir.resolve())
    write_config(config_path, cfg)

    print(f"✅ Project created at {project_dir}")
    return project_dir, config_path


if __name__ == "__main__":
    create_new_project()