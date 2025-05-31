import logging
from pathlib import Path

def setup_logger(name: str, log_file: str, log_dir: str = None, level=logging.INFO) -> logging.Logger:
    """
    Sets up a logger with the specified name and log file.

    This function creates a logger that writes log messages to both a file and the console.
    It ensures that duplicate handlers are not added to the logger. The logger is used to
    capture important events, errors, and debugging information during the execution of
    workflows in the bioinformatics pipeline.

    Args:
        name (str): Name of the logger. Typically, this is the name of the workflow or module.
        log_file (str): Name of the log file where messages will be stored.
        log_dir (str, optional): Directory where the log file will be saved. If not provided,
                                 the log file will be created in the current working directory.
        level (int): Logging level (e.g., logging.INFO, logging.DEBUG). Defaults to logging.INFO.

    Returns:
        logging.Logger: Configured logger instance.

    Usage in Pipelines:
        - Base Workflows:
            The base workflows (e.g., `BaseWorkflow`) typically use minimal logging, often
            relying on simple print statements for start and end messages. However, this
            logger can be integrated into base workflows to provide more structured logging
            if needed.

        - Non-Base Workflows:
            Advanced workflows (e.g., `CellRangerWorkflow`) use this logger extensively to
            capture detailed logs for each step of the workflow. This includes logging
            validation steps, configuration preparation, and execution of commands. The
            logger ensures that logs are written to both the console and a file, making it
            easier to debug and monitor the workflow.

    Example:
        ```python
        from utils.logger import setup_logger

        logger = setup_logger(
            name="CellRangerWorkflow",
            log_file="cell_ranger.log",
            log_dir="/path/to/logs",
            level=logging.DEBUG
        )
        logger.info("Workflow started.")
        logger.error("An error occurred.")
        ```
    """
    logger = logging.getLogger(name)

    # Avoid adding duplicate handlers
    if logger.hasHandlers():
        return logger

    logger.setLevel(level)

    # Determine the log file path
    if log_dir:
        log_dir_path = Path(log_dir)
        log_dir_path.mkdir(parents=True, exist_ok=True)  # Ensure the directory exists
        log_file_path = log_dir_path / log_file
    else:
        log_file_path = Path(log_file)

    # Create a file handler
    fh = logging.FileHandler(log_file_path)
    fh.setLevel(level)

    # Create a console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)

    # Create a formatter and set it for both handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    # Add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger

def get_log_file_path(base_dir: str) -> str:
    """
    Returns the path for the log file based on the base directory.

    This utility function is used to generate a consistent log file path for workflows.
    It is particularly useful for workflows that need to dynamically determine where
    to store their log files.

    Args:
        base_dir (str): The base directory where the log file should be stored.

    Returns:
        str: Full path to the log file.

    Usage in Pipelines:
        - Base Workflows:
            This function can be used to determine the log file path for base workflows,
            ensuring that logs are stored in a consistent location.

        - Non-Base Workflows:
            Advanced workflows can use this function to dynamically generate log file
            paths based on their specific requirements.

    Example:
        ```python
        from utils.logger import get_log_file_path

        base_dir = "/path/to/logs"
        log_file_path = get_log_file_path(base_dir)
        print(f"Log file will be stored at: {log_file_path}")
        ```
    """
    log_file_name = 'workflow_engine.log'
    return os.path.join(base_dir, log_file_name)