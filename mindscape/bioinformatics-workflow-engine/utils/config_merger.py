import yaml
from pathlib import Path
import tempfile

def merge_configs(default_config_path, project_config_path):
    """
    Merge default_config.yaml with project-specific config.yaml.

    Args:
        default_config_path (str or Path): Path to default_config.yaml.
        project_config_path (str or Path): Path to project-specific config.yaml.

    Returns:
        str: Path to the merged configuration file.
    """
    with open(default_config_path, 'r') as default_file:
        default_config = yaml.safe_load(default_file)

    with open(project_config_path, 'r') as project_file:
        project_config = yaml.safe_load(project_file)

    # Override relevant fields in default_config with values from project_config
    for key in ['project_name', 'experimenter', 'date', 'project_path', 'data_dir', 'results_dir', 'logs_dir', 'parameters']:
        if key in project_config:
            default_config[key] = project_config[key]

    # Write the merged configuration to a temporary file
    temp_config_file = tempfile.NamedTemporaryFile(delete=False, suffix=".yaml")
    with open(temp_config_file.name, 'w') as temp_file:
        yaml.dump(default_config, temp_file)

    return temp_config_file.name