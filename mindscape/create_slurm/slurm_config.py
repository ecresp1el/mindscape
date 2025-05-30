import yaml
from pathlib import Path

def collect_slurm_config():
    """
    Collect only the email for SLURM notifications. Use default parameters for everything else.

    Returns:
        dict: SLURM configuration.
    """
    print("🔧 Collecting SLURM configuration...")
    email = input("Enter your email for SLURM notifications (default: 'default@example.com'): ") or "default@example.com"

    # Default SLURM parameters
    slurm_config = {
        "account": "parent0",
        "time": "8:00:00",
        "memory": "32G",
        "cpus": 8,
        "mail-type": "FAIL",
        "mail-user": email,
    }
    print(f"Collected SLURM configuration: {slurm_config}")
    return slurm_config

def save_slurm_config(project_path, slurm_config):
    """
    Save SLURM configuration to a separate slurm_config.yaml file.

    Args:
        project_path (str or Path): Path to the project directory.
        slurm_config (dict): SLURM configuration.
    """
    slurm_config_path = Path(project_path) / "config/slurm_config.yaml"
    with open(slurm_config_path, "w") as file:
        yaml.safe_dump(slurm_config, file)

    print(f"✅ SLURM configuration saved to {slurm_config_path}")