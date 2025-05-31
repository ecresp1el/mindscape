from pathlib import Path
from mindscape.bioinformatics_workflow_engine.utils.slurm_job import SLURMJob

# ✅ POINT TO THE PROJECT PATH YOU WANT JOB FILES TO BE WRITTEN TO
project_path = Path("/nfs/turbo/umms-parent/Manny_test")

# ✅ Create a SLURM job (dry run)
job = SLURMJob(
    job_name="test_echo",
    command="echo 'Hello from SLURM dry run'",
    project_path=project_path,
    dry_run=True  # Only writes the .slurm script, doesn't submit
)

job.generate_script()