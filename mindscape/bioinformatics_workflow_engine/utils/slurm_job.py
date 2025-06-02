from pathlib import Path
import uuid
import re

class SLURMJob:
    def __init__(
        self,
        job_name: str,
        command: str,
        project_path: Path,
        time: str = "12:00:00",
        mem: str = "64G",
        cpus: int = 8,
        account: str = "parent0",
        email: str = "elcrespo@umich.edu",
        dry_run: bool = True,
    ):
        self.job_name = job_name
        self.command = command
        self.project_path = Path(project_path)
        self.time = time
        self.mem = mem
        self.cpus = cpus
        self.account = account
        self.email = email
        self.dry_run = dry_run
        self.time = self._validate_time_format(str(self.time))
        self.slurm_dir = self.project_path / "logs" / "slurm_scripts"
        self.log_dir = self.project_path / "logs" / "slurm_logs"
        self.script_path = self.slurm_dir / f"{job_name}_{uuid.uuid4().hex[:6]}.slurm" # Unique script name to avoid collisions

    def generate_script(self):
        self.slurm_dir.mkdir(parents=True, exist_ok=True)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        script = f"""#!/bin/bash

# ----------------------------------------------------------------------------
# SLURM Configuration
# ----------------------------------------------------------------------------
#SBATCH --job-name={self.job_name}
#SBATCH --output={self.log_dir}/{self.job_name}_{self.script_path.stem}.out
#SBATCH --error={self.log_dir}/{self.job_name}_{self.script_path.stem}.err
#SBATCH --account={self.account}
#SBATCH --time={self.time}
#SBATCH --mem={self.mem}
#SBATCH --cpus-per-task={self.cpus}
#SBATCH --mail-type=FAIL
#SBATCH --mail-user={self.email}

# ----------------------------------------------------------------------------
# Environment and software setup
# ----------------------------------------------------------------------------
source ~/.bashrc
conda activate mindscape-env

# ----------------------------------------------------------------------------
# Command to run
# ----------------------------------------------------------------------------
set -e
{self.command}
"""
        self.script_path.write_text(script)
        print(f"[Dry Run] SLURM script written to: {self.script_path}")
        print("----------- SLURM SCRIPT PREVIEW -----------")
        print(script)
        print("--------------------------------------------")
        return self.script_path

    def submit(self):
        if self.dry_run:
            self.generate_script()
            print(f"[Dry Run] SLURM script generated for {self.job_name} (not submitted)")
            return

        self.generate_script()
        import subprocess
        result = subprocess.run(["sbatch", str(self.script_path)], capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"SLURM submission failed:\n{result.stderr}")
        job_id = result.stdout.strip().split()[-1]
        print(f"✅ SLURM job submitted: {self.job_name} (Job ID: {job_id})")
        return job_id

    def _validate_time_format(self, time_str):
        """
        Validates and auto-corrects SLURM time format.
        Accepts formats like 'HH:MM:SS' or integers in hours (e.g., '12').
        """
        if re.fullmatch(r"\d{1,2}:\d{2}:\d{2}", time_str):
            return time_str
        elif time_str.isdigit():
            return f"{time_str.zfill(2)}:00:00"
        else:
            raise ValueError(f"Invalid SLURM time format: {time_str}. Use 'HH:MM:SS' or an integer representing hours.")