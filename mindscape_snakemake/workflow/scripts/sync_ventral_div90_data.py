# workflow/scripts/sync_ventral_div90_data.py
import os
import shutil
import sys
import yaml
from pathlib import Path

def sync_data(config_path):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    project_path = Path(config["project_path"])
    source_dir = Path("/nfs/turbo/umms-parent/Manny_test/10496-MW-reanalysis/outs/per_sample_outs")
    dest_dir = project_path / "data" / "ventral_sosrs_div90"
    dest_dir.mkdir(parents=True, exist_ok=True)

    sample_dirs = sorted([p for p in source_dir.iterdir() if p.is_dir()])

    for sample_dir in sample_dirs:
        sample_id = sample_dir.name
        count_path = sample_dir / "count"
        if count_path.exists():
            target_path = dest_dir / sample_id
            print(f"🔁 Copying {count_path} → {target_path}")
            shutil.copytree(count_path, target_path, dirs_exist_ok=True)
        else:
            print(f"⚠️ Skipping {sample_id}: count/ folder not found.")

    print("✅ Sync complete.")

# Support both standalone and Snakemake usage
if __name__ == "__main__":
    if len(sys.argv) > 1:
        sync_data(sys.argv[1])
    else:
        from snakemake import snakemake  # noqa: F401
        sync_data(snakemake.input[0])