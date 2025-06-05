# workflow/scripts/sync_ventral_div90_data.py
import shutil
from pathlib import Path
import yaml

# This script uses Snakemake's injected `snakemake` object

# Load config to resolve paths
with open(snakemake.input[0]) as f:
    config = yaml.safe_load(f)

source_root = Path(config["shared_turbo_directory"]) / "Manny_test" / "10496-MW-reanalysis" / "outs" / "per_sample_outs"
dest_base = Path(config["project_path"]) / "data" / "ventral_sosrs_div90"

# Ensure destination exists
dest_base.mkdir(parents=True, exist_ok=True)

# For each subdirectory (i.e., sample_id)
for sample_dir in source_root.iterdir():
    if sample_dir.is_dir():
        sample_id = sample_dir.name
        src = sample_dir / "count"
        dest = dest_base / sample_id
        if src.exists():
            print(f"📁 Copying from {src} to {dest}")
            shutil.copytree(src, dest, dirs_exist_ok=True)

# Create marker file to satisfy Snakemake
touch = Path(snakemake.output[0])
touch.parent.mkdir(parents=True, exist_ok=True)
touch.touch()
