def submit_slurm_job(script_path: str, job_name: str, partition: str = "default", time: str = "01:00:00", nodes: int = 1, cpus: int = 1) -> str:
    """
    Submits a job to SLURM.

    Parameters
    ----------
    script_path : str
        Path to the script to be executed.
    job_name : str
        Name of the job.
    partition : str, optional
        SLURM partition to submit the job to (default is "default").
    time : str, optional
        Time limit for the job (default is "01:00:00").
    nodes : int, optional
        Number of nodes to allocate (default is 1).
    cpus : int, optional
        Number of CPUs to allocate (default is 1).

    Returns
    -------
    str
        Job ID of the submitted job.
    """
    command = f"sbatch --job-name={job_name} --partition={partition} --time={time} --nodes={nodes} --cpus-per-task={cpus} {script_path}"
    job_id = os.popen(command).read().strip()
    return job_id


def check_job_status(job_id: str) -> str:
    """
    Checks the status of a SLURM job.

    Parameters
    ----------
    job_id : str
        The job ID to check.

    Returns
    -------
    str
        The status of the job.
    """
    command = f"squeue --job {job_id} --noheader"
    status = os.popen(command).read().strip()
    return status if status else "Job not found or completed."


def cancel_slurm_job(job_id: str) -> None:
    """
    Cancels a SLURM job.

    Parameters
    ----------
    job_id : str
        The job ID to cancel.
    """
    command = f"scancel {job_id}"
    os.system(command)