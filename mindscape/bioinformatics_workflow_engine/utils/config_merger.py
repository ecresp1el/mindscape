import yaml
from pathlib import Path

def merge_configs(default_config_path, project_config_path, output_path=None):
    """
    Merge default_config.yaml with project-specific config.yaml.

    Args:
        default_config_path (str or Path): Path to default_config.yaml.
        project_config_path (str or Path): Path to project-specific config.yaml.
        output_path (str or Path, optional): Path to save the merged configuration. Defaults to project config directory.

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

    # Determine the output path
    if output_path is None:
        project_dir = Path(project_config_path).parent
        output_path = project_dir / "merged_config.yaml"

    # Save the merged configuration to the output path
    with open(output_path, 'w') as output_file:
        yaml.dump(default_config, output_file)

    return str(output_path)