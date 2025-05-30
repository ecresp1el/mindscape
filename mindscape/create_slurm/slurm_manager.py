import subprocess
import logging
from pathlib import Path

class SlurmManager:
    """
    Handles SLURM job submission and resource allocation.
    """

    def __init__(self, slurm_config, log_dir=None):
        """
        Initialize the SlurmManager with SLURM-specific configurations.

        Args:
            slurm_config (dict): SLURM configuration (e.g., memory, cores, partition).
            log_dir (str or Path, optional): Directory for SLURM log files. Defaults to None.
        """
        self.slurm_config = slurm_config
        self.log_dir = Path(log_dir) if log_dir else None
        self.logger = logging.getLogger("SlurmManager")

        # Ensure the log directory exists
        if self.log_dir:
            self.log_dir.mkdir(parents=True, exist_ok=True)

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
            f"--cpus-per-task={self.slurm_config.get('cpus', 1)}",
            f"--mem={self.slurm_config.get('memory', '4G')}",
            f"--partition={self.slurm_config.get('partition', 'default')}",
            f"--time={self.slurm_config.get('time', '1:00:00')}",
            f"--mail-type={self.slurm_config.get('mail-type', 'FAIL')}",
            f"--mail-user={self.slurm_config.get('mail-user', 'default@example.com')}",
        ]

        # Add log file paths
        if self.log_dir:
            slurm_command.append(f"--output={self.log_dir}/{job_name}_{pipeline_step}.out")
            slurm_command.append(f"--error={self.log_dir}/{job_name}_{pipeline_step}.err")

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