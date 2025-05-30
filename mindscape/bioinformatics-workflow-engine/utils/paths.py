from pathlib import Path

def get_project_root() -> Path:
    """Returns the root directory of the project."""
    return Path(__file__).resolve().parent.parent

def get_data_dir() -> Path:
    """Returns the path to the data directory."""
    return get_project_root() / 'data'

def get_results_dir() -> Path:
    """Returns the path to the results directory."""
    return get_project_root() / 'results'

def get_logs_dir() -> Path:
    """Returns the path to the logs directory."""
    return get_project_root() / 'logs'

def get_config_dir() -> Path:
    """Returns the path to the config directory."""
    return get_project_root() / 'config'

def get_workflows_dir() -> Path:
    """Returns the path to the pipelines directory."""
    return get_project_root() / 'pipelines'