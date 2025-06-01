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

    This class is designed to serve as a "blueprint" for creating workflows in the bioinformatics workflow engine.
    It provides common functionality that all workflows can use, such as loading configuration files, setting up
    project paths, and logging the start and end of a workflow.

    Think of this class as a template. Other workflows (like CellRangerWorkflow or QCWorkflow) will "inherit" from
    this class, meaning they will reuse the functionality provided here while adding their own specific behavior.

    Key Features:
    - **Configuration Management**: Automatically loads a configuration file (YAML format) for the workflow.
    - **Path Setup**: Ensures that required directories for the workflow exist.
    - **Logging**: Provides simple methods to log when a workflow starts and ends.
    - **Extensibility**: Requires subclasses to define their own `run` method, which contains the specific steps
      for that workflow.

    Why Use This Class?
    - It avoids repeating the same code in every workflow.
    - It ensures all workflows follow the same structure, making the system easier to maintain and extend.

    Example:
        class MyWorkflow(BaseWorkflow):
            def run(self):
                self.log_start()
                # Perform workflow-specific tasks here
                self.log_end()
    """

    def __init__(self, config_path):
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
        self.config_path = Path(config_path).resolve()  # Convert the config path to an absolute path
        self.config = self.load_config()  # Load the configuration file
        self.workflow_name = self.__class__.__name__  # Get the name of the workflow (e.g., "CellRangerWorkflow")
        # This line dynamically retrieves the name of the class to which the current object belongs.
        # It uses the `__class__` attribute of the object to access its class and then retrieves the
        # `__name__` attribute of the class, which contains the name of the class as a string.
        # For example, if the current object is an instance of the `CellRangerWorkflow` class,
        # this line will set `self.workflow_name` to the string "CellRangerWorkflow".
        #
        # Why is this useful?
        # - It allows the workflow to "know" its own name, which can be used for logging, debugging,
        #   or dynamically identifying the type of workflow at runtime.
        # - This approach ensures that the workflow name is always consistent with the class name,
        #   even if the class is renamed or subclassed in the future.
        # - It avoids hardcoding the workflow name, making the code more flexible and maintainable.
        
        # SLRUM resources configuration, pulls cpu, memory, and time from the config file
        # defaults are set if not specified in the config file 
        slurm_cfg = self.config.get("slurm", {}) # Default to an empty dictionary if 'slurm' is not defined
        self.slurm_cpus = slurm_cfg.get("cpus", 8) # Default to 8 CPUs if not specified
        self.slurm_mem = slurm_cfg.get("mem", "32G") # Default to 32G memory if not specified
        self.slurm_time = slurm_cfg.get("time", "08:00:00") # Default to 8 hours if not specified

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
        self.project_path = Path(self.config.get('project_path', '/tmp/default_project'))
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
        print(f"Starting workflow: {self.workflow_name}")

    def log_end(self):
        """
        Log the end of the workflow.

        This method prints a message to indicate that the workflow has completed. It is useful for tracking
        the progress of workflows, especially when running multiple workflows in a pipeline.

        Example Output:
            If the workflow is named "QCWorkflow", this method will print:
            "Completed workflow: QCWorkflow"
        """
        print(f"Completed workflow: {self.workflow_name}")

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
            
            print(f"[Dry Run] Generating SLURM job: {job_name} with the user email {email}")
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
            print(f"[Dry Run] Would run locally: {command}")
            if not dry_run:
                subprocess.run(command, shell=True, check=True)
            return None