import logging
import os
from pathlib import Path

def setup_logger(name: str, log_file: str, log_dir: str = None, level=logging.INFO) -> logging.Logger:
    """
    Sets up a logger with the specified name and log file.

    Args:
        name (str): Name of the logger.
        log_file (str): Name of the log file.
        log_dir (str, optional): Directory where the log file will be saved. Defaults to None.
        level (int): Logging level. Defaults to logging.INFO.

    Returns:
        logging.Logger: Configured logger.
    """
    logger = logging.getLogger(name)
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
    if not logger.handlers:  # Avoid adding duplicate handlers
        logger.addHandler(fh)
        logger.addHandler(ch)

    return logger

def get_log_file_path(base_dir: str) -> str:
    """Returns the path for the log file based on the base directory."""
    log_file_name = 'workflow_engine.log'
    return os.path.join(base_dir, log_file_name)