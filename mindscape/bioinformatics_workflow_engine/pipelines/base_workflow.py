from __future__ import annotations
import os
import logging
from pathlib import Path
import yaml
from mindscape.bioinformatics_workflow_engine.utils.slurm_job import SLURMJob
import subprocess

class BaseWorkflow:
    """
    BaseWorkflow: A Foundation for All Workflows
    
    BaseWorkflow is a superclass that provides essential functionality for all workflows in the bioinformatics workflow engine.

    This class is designed to serve as a "blueprint" for creating workflows in the bioinformatics workflow engine.
    It provides common functionality that all workflows can use, such as loading configuration files, setting up
    project paths, and logging the start and end of a workflow.

    Workflow Completion Tracking:
    - Workflows that complete successfully mark themselves using a `.completed` file located in:
        <project_path>/logs/<workflow_name>.completed
    - When `run_workflows.py` is executed, it uses `is_already_completed()` to skip rerunning steps
      that have already finished.
    - If you need to force re-run, delete the `.completed` file manually or implement a `force_rerun` flag check.

    What Subclasses Inherit:
    ------------------------
    Subclasses of BaseWorkflow inherit the following behavior and attributes:

    - Initialization:
        - Automatically loads the YAML configuration file specified by `config_path`.
        - Sets up a scoped logger named after the subclass (e.g., "workflow.MyWorkflow").
        - Parses SLURM resource settings (`cpus`, `mem`, `time`) from the configuration with sensible defaults.
        - Sets `self.workflow_name` to the subclass's class name.
    
    - Attributes available to subclasses:
        - `self.config`: The loaded configuration dictionary.
        - `self.project_path`: The root project path as specified in the config (after calling `setup_paths()`).
        - `self.logger`: A logger instance scoped to the workflow subclass.
        - SLURM-related attributes: `self.slurm_cpus`, `self.slurm_mem`, `self.slurm_time`.
    
    - Methods provided:
        - Logging helpers: `log_start()`, `log_end()`, and `mark_in_progress()` for consistent workflow status logging.
        - Status tracking: `is_already_completed()`, `mark_completed()`, `mark_failed(reason)`, and `get_status()`.
        - Path management: `setup_paths()` to create necessary directories.
        - Job submission: `submit_job(command, job_name, dry_run=True)` to submit commands either locally or via SLURM.
        - Configuration loading: `load_config()` for reading YAML config files.
    
    - Abstract method:
        - `run()`: Must be implemented by subclasses to define the specific workflow steps.
    
    Subclasses should call `setup_paths()` early in their workflow to ensure directories exist, use the logging methods
    to track progress, and implement their own `run()` method to execute workflow-specific logic.
    
    Example:
        class MyWorkflow(BaseWorkflow):
            def run(self):
                self.log_start()
                self.setup_paths()
                # Workflow-specific tasks here
                self.log_end()
    """

    def __init__(self, config_path, logger=None):
        """
        Initialize the BaseWorkflow.

        When a workflow is created, it needs to know where its configuration file is located.
        This method sets up the workflow by:
        1. Resolving the configuration file path to an absolute path.
        2. Loading the configuration file into memory.
        3. Storing the name of the workflow (based on the class name).

        Args:
            config_path (str): The path to the YAML configuration file for the workflow.
        """
        # Assign a logger scoped specifically to the subclass of the workflow.
        # This allows each workflow to have its own logging context, making logs easier to filter and understand.
        # If a logger is provided externally, use it; otherwise, create one using the class name.
        self.logger = logger or logging.getLogger(f"workflow.{self.__class__.__name__}")

        # Convert the provided config_path to an absolute Path object.
        # This ensures consistent access to the configuration file regardless of the current working directory.
        self.config_path = Path(config_path).resolve()

        # Load the YAML configuration from the config_path into a dictionary.
        # This config will be used throughout the workflow to access settings and parameters.
        self.config = self.load_config()

        # Store the name of the workflow as the class name of the subclass.
        # This supports consistent naming across logs, file markers, and job submissions,
        # and avoids hardcoding the workflow name in multiple places.
        self.workflow_name = self.__class__.__name__

        # Extract SLURM-related resource parameters from the configuration.
        # Provide sensible defaults if these parameters are not specified:
        # - cpus: number of CPU cores to request (default: 8)
        # - mem: amount of memory to allocate (default: "32G")
        # - time: maximum runtime for the job (default: "08:00:00")
        slurm_cfg = self.config.get("slurm", {})
        self.slurm_cpus = int(slurm_cfg.get("cpus", 8))  # Ensure CPUs is an integer
        self.slurm_mem = str(slurm_cfg.get("mem", "32G"))  # Ensure memory is a string
        self.logger.debug(f"Parsed SLURM resources: {self.slurm_cpus} CPUs, {self.slurm_mem}")

        raw_time = slurm_cfg.get("time", "08:00:00")
        self.logger.debug(f"Raw SLURM time: {raw_time} (type: {type(raw_time)})")

        # Normalize the time format to HH:MM:SS if an integer or digit string is provided.
        # This allows users to specify time as an integer number of hours for convenience.
        if isinstance(raw_time, int):
            self.slurm_time = f"{str(raw_time).zfill(2)}:00:00"
        elif isinstance(raw_time, str) and raw_time.isdigit():
            self.slurm_time = f"{raw_time.zfill(2)}:00:00"
        else:
            self.slurm_time = raw_time

        self.logger.debug(f"Parsed SLURM time: {self.slurm_time} (type: {type(self.slurm_time)})")

    def load_config(self):
        """
        Load the YAML configuration file.

        This method reads the configuration file specified during initialization and converts it into a Python
        dictionary. The configuration contains important settings that the workflow needs to run.

        Returns:
            dict: The contents of the configuration file as a dictionary.

        Example:
            If the YAML file contains:
            ```
            project_path: /home/user/project
            workflows:
              - name: QCWorkflow
                enabled: true
            ```
            Then this method will return:
            {
                "project_path": "/home/user/project",
                "workflows": [
                    {"name": "QCWorkflow", "enabled": True}
                ]
            }
        """
        with open(self.config_path, 'r') as file:
            return yaml.safe_load(file)

    def setup_paths(self):
        """
        Set up the required project paths.

        Many workflows need to save results or logs in specific directories. This method ensures that the
        required directories exist. If they don't, it creates them.

        Example:
            If the configuration specifies:
            ```
            project_path: /home/user/project
            ```
            This method will ensure that `/home/user/project` exists. If it doesn't, it will create it.

        Why This is Useful:
        - It prevents errors caused by missing directories.
        - It ensures all workflows follow the same directory structure.
        """
        if "project_path" not in self.config:
            raise KeyError("❌ 'project_path' is missing from the configuration file.")
        self.project_path = Path(self.config["project_path"])
        if not self.project_path.exists():
            self.project_path.mkdir(parents=True, exist_ok=True)

    def log_start(self):
        """
        Log the start of the workflow.

        This method prints a message to indicate that the workflow has started. It is useful for tracking
        the progress of workflows, especially when running multiple workflows in a pipeline.

        Example Output:
            If the workflow is named "QCWorkflow", this method will print:
            "Starting workflow: QCWorkflow"
        """
        self.logger.info(f"🟢 Starting workflow: {self.workflow_name}")

    def mark_in_progress(self):
        """
        Marks this workflow as in progress. Creates a .in_progress file.
        This file is NOT removed automatically when completed — it serves as a record that the workflow began.
        """
        marker = self.get_completion_marker_path().with_suffix(".in_progress")
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text("IN PROGRESS\n")
        self.logger.info(f"📍 Marked workflow as in progress: {self.workflow_name}")
        
        
    def get_completion_marker_path(self) -> Path:
        """
        Return the path to the marker file used to indicate this workflow is completed.
        """
        # The marker file is stored under <project_path>/logs/<workflow_name>.completed
        return self.project_path / "logs" / f"{self.workflow_name}.completed"

    def is_already_completed(self) -> bool:
        """
        Check whether this workflow has already been completed by looking for a .completed marker file.
        If 'force_rerun' is True in config, always return False to force re-execution.
        """
        if self.config.get("force_rerun", False):
            return False
        return self.get_completion_marker_path().exists()

    def mark_completed(self):
        """
        Marks this workflow as completed. Leaves the .in_progress file intact for auditability.
        """
        completed = self.get_completion_marker_path()
        completed.parent.mkdir(parents=True, exist_ok=True)
        completed.write_text("COMPLETED\n")
        self.logger.info(f"🏁 Marked workflow as completed: {self.workflow_name}")

    def mark_failed(self, reason: str = "Unspecified"):
        """
        Marks this workflow as failed by writing a .failed file with the reason.
        """
        failed_marker = self.get_completion_marker_path().with_suffix(".failed")
        failed_marker.parent.mkdir(parents=True, exist_ok=True)
        failed_marker.write_text(f"FAILED: {reason}\n")
        self.logger.error(f"❌ Marked workflow as failed: {self.workflow_name} — {reason}")

    def get_status(self) -> str:
        """
        Returns the current status of the workflow based on marker files.
        """
        if self.get_completion_marker_path().exists():
            return "completed"
        elif self.get_completion_marker_path().with_suffix(".in_progress").exists():
            return "in_progress"
        elif self.get_completion_marker_path().with_suffix(".failed").exists():
            return "failed"
        return "not_started"

    def log_end(self):
        """
        Log the end of the workflow and always mark it as completed.
        """
        self.logger.info(f"✅ Completed workflow: {self.workflow_name}")
        self.mark_completed()

    def run(self):
        """
        Run the workflow.

        This method is a placeholder that must be implemented by subclasses. Each specific workflow
        (e.g., CellRangerWorkflow, QCWorkflow) will define its own method, which contains the steps
        required to execute that workflow.

        Why This is Important:
        - It ensures that every workflow defines its own behavior.
        - It provides a consistent interface for running workflows, making it easier to manage them.

        Example:
            class QCWorkflow(BaseWorkflow):
                def run(self):
                    self.log_start()
                    # Perform quality control tasks
                    self.log_end()

        Raises:
            NotImplementedError: If a subclass does not implement this method.
        """
        raise NotImplementedError("Subclasses must implement the run method.")
    
    def submit_job(self, command: str, job_name: str, dry_run: bool = True):
        """
        Submit a command either via SLURM or locally, based on self.use_slurm and dry_run.
        """
        if getattr(self, "use_slurm", False):
            
            email = self.config.get("email", "elcrespo@umich.edu") # Default email if not specified
            self.logger.info(f"[Dry Run] Generating SLURM job: {job_name} with email: {email}")
            # Create a SLURM job instance
            job = SLURMJob(
                job_name=job_name,
                command=command,
                email=email,
                project_path=self.project_path,
                dry_run=dry_run, 
                cpus=self.slurm_cpus,
                mem=self.slurm_mem, #i.e. 64G so remember to use the string format within 
                time=self.slurm_time
            )
            return job.generate_script() if dry_run else job.submit()
        else:
            self.logger.info(f"[Dry Run] Would run locally: {command}")
            if not dry_run:
                self.logger.info(f"Running command locally: {command}")
                subprocess.run(command, shell=True, check=True)
            return None