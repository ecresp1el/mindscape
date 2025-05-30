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
        slurm_config_path = Path(project_path) / "config/slurm_config.yaml"
        with open(slurm_config_path, "r") as file:
            self.slurm_config = yaml.safe_load(file)

        self.log_dir = Path(log_dir) if log_dir else Path(project_path) / "logs"
        self.log_dir.mkdir(parents=True, exist_ok=True)

        self.logger = logging.getLogger("SlurmManager")

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
        slurm_command = [
            "sbatch",
            f"--job-name={job_name}_{pipeline_step}",
            f"--account={self.slurm_config.get('account', 'parent0')}",
            f"--time={self.slurm_config.get('time', '8:00:00')}",
            f"--mem={self.slurm_config.get('memory', '32G')}",
            f"--cpus-per-task={self.slurm_config.get('cpus', 8)}",
            f"--mail-type={self.slurm_config.get('mail-type', 'FAIL')}",
            f"--mail-user={self.slurm_config.get('mail-user', 'default@example.com')}",
            f"--output={self.log_dir}/{job_name}_{pipeline_step}.out",
            f"--error={self.log_dir}/{job_name}_{pipeline_step}.err",
        ]

        # Add the command to execute
        slurm_command.append("--wrap")
        slurm_command.append(command)

        self.logger.info(f"Submitting SLURM job: {' '.join(slurm_command)}")
        result = subprocess.run(slurm_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if result.returncode != 0:
            self.logger.error(f"SLURM job submission failed: {result.stderr}")
            raise RuntimeError(f"SLURM job submission failed: {result.stderr}")

        job_id = result.stdout.strip().split()[-1]
        self.logger.info(f"SLURM job submitted successfully. Job ID: {job_id}")
        return job_id