from snakemake.io import glob_wildcards
import os

def get_input_files(input_pattern: str) -> list:
    """
    Returns a list of input files matching the given pattern.
    
    Parameters
    ----------
    input_pattern : str
        The glob pattern to match input files.

    Returns
    -------
    list
        A list of matched input file paths.
    """
    return glob_wildcards(input_pattern).input

def create_snakemake_rule(rule_name: str, input_files: list, output_files: list, shell_command: str):
    """
    Creates a Snakemake rule.

    Parameters
    ----------
    rule_name : str
        The name of the Snakemake rule.
    input_files : list
        List of input files for the rule.
    output_files : list
        List of output files for the rule.
    shell_command : str
        The shell command to execute for the rule.
    """
    rule rule_name:
        input: input_files
        output: output_files
        shell:
            shell_command

def run_snakemake(workflow_file: str):
    """
    Executes the Snakemake workflow defined in the specified file.

    Parameters
    ----------
    workflow_file : str
        The path to the Snakemake workflow file.
    """
    os.system(f"snakemake -s {workflow_file}")