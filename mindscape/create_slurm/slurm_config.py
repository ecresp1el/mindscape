import yaml
from pathlib import Path
import sys

def collect_slurm_config():
    print("DEBUG: Collecting SLURM configuration...")
    
    if sys.stdin.isatty():
        email = input("Enter your email for SLURM notifications (default: 'default@example.com'): ") or "default@example.com"
    else:
        print("⚠️ No interactive terminal detected. Using default email.")
        email = "default@example.com"

    slurm_config = {
        "account": "parent0",
        "time": "8:00:00",
        "memory": "32G",
        "cpus": 8,
        "mail-type": "FAIL",
        "mail-user": email,
    }
    print(f"DEBUG: Collected SLURM configuration: {slurm_config}")
    return slurm_config

def save_slurm_config(project_path, slurm_config):
    print(f"DEBUG: Saving SLURM configuration to {project_path}/config/slurm_config.yaml...")
    slurm_config_path = Path(project_path) / "config/slurm_config.yaml"
    with open(slurm_config_path, "w") as file:
        yaml.safe_dump(slurm_config, file)

    print(f"DEBUG: SLURM configuration saved to {slurm_config_path}")