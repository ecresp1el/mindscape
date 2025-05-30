import yaml

def collect_slurm_config():
    """
    Collect SLURM-specific details from the user.

    Returns:
        dict: SLURM configuration.
    """
    print("🔧 Collecting SLURM configuration...")
    account = input("Enter your SLURM account (default: 'default_account'): ") or "default_account"
    email = input("Enter your email for SLURM notifications (default: 'default@example.com'): ") or "default@example.com"
    partition = input("Enter the SLURM partition (default: 'default'): ") or "default"
    time = input("Enter the default time allocation (e.g., '1:00:00', default: '1:00:00'): ") or "1:00:00"
    memory = input("Enter the default memory allocation (e.g., '4G', default: '4G'): ") or "4G"
    cpus = input("Enter the default number of CPUs (default: 1): ") or 1

    slurm_config = {
        "account": account,
        "mail-user": email,
        "partition": partition,
        "time": time,
        "memory": memory,
        "cpus": int(cpus),
        "mail-type": "FAIL",
    }
    print(f"Collected SLURM configuration: {slurm_config}")
    return slurm_config

def save_slurm_config(config_path, slurm_config):
    print(f"Saving SLURM configuration to {config_path}...")
    with open(config_path, "r") as file:
        config = yaml.safe_load(file)

    config["slurm"] = slurm_config

    with open(config_path, "w") as file:
        yaml.safe_dump(config, file)

    print(f"✅ SLURM configuration saved: {slurm_config}")