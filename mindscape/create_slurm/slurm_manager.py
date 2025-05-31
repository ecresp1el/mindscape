import yaml
from pathlib import Path
import subprocess
import logging

class SlurmManager:
    """
    Handles SLURM job submission and resource allocation.
    """

    def __init__(self, project_path, log_dir=None):
        """
        Initialize the SlurmManager with SLURM-specific configurations.

        Args:
            project_path (str or Path): Path to the project directory.
            log_dir (str or Path, optional): Directory for SLURM log files. Defaults to None.
        """
        self.project_path = Path(project_path)
        self.slurm_config_path = self.project_path / "config/slurm_config.yaml"
        self.log_dir = Path(log_dir) if log_dir else self.project_path / "logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self.logger = logging.getLogger("SlurmManager")
        self.load_slurm_config()

    def load_slurm_config(self):
        """Load the SLURM configuration from the project folder."""
        if not self.slurm_config_path.exists():
            raise FileNotFoundError(f"SLURM configuration file not found at {self.slurm_config_path}")
        with open(self.slurm_config_path, "r") as file:
            self.slurm_config = yaml.safe_load(file)
        self.logger.info(f"Loaded SLURM configuration: {self.slurm_config}")

    def submit_job(self, command, job_name, pipeline_step):
        """
        Submit a job to SLURM.

        Args:
            command (str): The command to execute.
            job_name (str): The name of the SLURM job.
            pipeline_step (str): The pipeline substep name.

        Returns:
            str: SLURM job ID.
        """
        print("🔍 [DEBUG] Loading SLURM configuration...")
        self.load_slurm_config()
        print(f"✅ [DEBUG] SLURM configuration loaded: {self.slurm_config}")

        setup_env = """
        # ------------------------------------------------------------------------------
        # Environment and software setup
        # ------------------------------------------------------------------------------

        # Initialize Conda (needed to use 'conda activate')
        source ~/.bashrc

        # Load required UMich modules first (do not override Conda Python later)
        module purge
        module load Bioinformatics cellranger

        # Activate your fully configured MindScape Conda environment
        conda activate mindscape-env

        # Ensure Python can find the MindScape package
        export PYTHONPATH=/home/elcrespo/Desktop/githubprojects/MindScape

        # Sanity checks
        echo "✅ Python path: $(which python)"
        python --version
        echo "✅ Cell Ranger path: $(which cellranger)"
        cellranger --version
        """

        print("🔍 [DEBUG] Preparing full command for SLURM...")
        full_command = f"{setup_env}\n{command}"
        print(f"✅ [DEBUG] Full command prepared:\n{full_command}")

        # Write the SLURM script to a file
        slurm_script_path = self.log_dir / f"{job_name}_{pipeline_step}.slurm"
        with open(slurm_script_path, "w") as script_file:
            script_file.write("#!/bin/bash\n")
            script_file.write(f"#SBATCH --job-name={job_name}_{pipeline_step}\n")
            script_file.write(f"#SBATCH --account={self.slurm_config.get('account', 'parent0')}\n")
            script_file.write(f"#SBATCH --time={self.slurm_config.get('time', '8:00:00')}\n")
            script_file.write(f"#SBATCH --mem={self.slurm_config.get('memory', '32G')}\n")
            script_file.write(f"#SBATCH --cpus-per-task={self.slurm_config.get('cpus', 8)}\n")
            script_file.write(f"#SBATCH --mail-type={self.slurm_config.get('mail-type', 'FAIL')}\n")
            script_file.write(f"#SBATCH --mail-user={self.slurm_config.get('mail-user', 'default@example.com')}\n")
            script_file.write(f"#SBATCH --output={self.log_dir}/{job_name}_{pipeline_step}.out\n")
            script_file.write(f"#SBATCH --error={self.log_dir}/{job_name}_{pipeline_step}.err\n")
            script_file.write("\n")
            script_file.write(setup_env)
            script_file.write("\n")
            script_file.write(command)
            script_file.write("\n")

        # Submit the script via sbatch
        print("🔍 [DEBUG] Submitting SLURM job script...")
        result = subprocess.run(["sbatch", str(slurm_script_path)],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        print("🔍 [DEBUG] SLURM submission result:")
        print(f"STDOUT:\n{result.stdout}")
        print(f"STDERR:\n{result.stderr}")

        if result.returncode != 0:
            print("❌ [ERROR] SLURM job submission failed!")
            raise RuntimeError(f"SLURM job submission failed: {result.stderr}")

        job_id = result.stdout.strip().split()[-1]


        return job_id