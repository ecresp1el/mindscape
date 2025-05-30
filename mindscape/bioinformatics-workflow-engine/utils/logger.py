import logging
import os

def setup_logger(name: str, log_file: str, level=logging.INFO) -> logging.Logger:
    """Sets up a logger with the specified name and log file."""
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Create a file handler
    fh = logging.FileHandler(log_file)
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
    """Returns the path for the log file based on the base directory."""
    log_file_name = 'workflow_engine.log'
    return os.path.join(base_dir, log_file_name)