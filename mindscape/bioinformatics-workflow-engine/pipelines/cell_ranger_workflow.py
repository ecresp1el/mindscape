import yaml
from pathlib import Path
import shutil
import subprocess
import tempfile
from .base_workflow import BaseWorkflow
from mindscape.create_slurm.slurm_manager import SlurmManager

class CellRangerWorkflow(BaseWorkflow):
    def __init__(self, config):
        super().__init__(config)

        # Load the configuration if a path is provided
        if isinstance(config, str) or isinstance(config, Path):
            with open(config, "r") as file:
                config = yaml.safe_load(file)

        self.workflow_name = "CellRangerWorkflow"

        # Base paths and subpaths
        self.project_path = Path(config["project_path"])
        self.results_dir = self.project_path / "results" / f"cellranger_multi_{self.workflow_name}"
        self.logs_dir = self.project_path / "logs"

        # Reference genome and other files
        self.turbo_ref_base = "/nfs/turbo/umms-parent/10X_Human_Refs"
        self.ref_subpath = "2020-A/Ref_genome/refdata-gex-GRCh38-2020-A"
        self.ref_genome = Path(self.turbo_ref_base) / self.ref_subpath

        self.multi_config_source = "/nfs/turbo/umms-parent/Accessible_multi-config_csvs/Carmen_Miranda_scRNAseq /90 Day results/fastqs_10496-MW/multi-config.csv"
        self.probe_file = "/nfs/turbo/umms-parent/10X_Human_Refs/2020-A/Probe_set/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv"
        self.output_id = "10496-MW-reanalysis"

    def validate_paths(self):
        """Validate that required paths exist."""
        # Validate reference genome
        if not self.ref_genome.exists():
            raise FileNotFoundError(f"❌ ERROR: Reference folder not found at {self.ref_genome}")
        print(f"🧬 Using reference genome: {self.ref_genome}")

        # Validate multi_config.csv source
        print(f"DEBUG: Checking multi-config source path: {self.multi_config_source}")
        if not Path(self.multi_config_source).exists():
            raise FileNotFoundError(f"❌ ERROR: Multi-config file not found at {self.multi_config_source}")
        print(f"📄 Multi-config source found: {self.multi_config_source}")

        # Validate probe file
        if not Path(self.probe_file).exists():
            raise FileNotFoundError(f"❌ ERROR: Probe file not found at {self.probe_file}")
        print(f"🧬 Probe file found: {self.probe_file}")

    def prepare_multi_config(self):
        """Prepare the multi_config.csv file."""
        config_dest = self.results_dir / "multi_config.csv"

        # Ensure the results directory exists
        print(f"Ensuring results directory exists: {self.results_dir}")
        self.results_dir.mkdir(parents=True, exist_ok=True)

        # Copy the multi_config.csv file
        print(f"Copying multi_config.csv from {self.multi_config_source} to {config_dest}")
        try:
            shutil.copy(self.multi_config_source, config_dest)
        except Exception as e:
            raise RuntimeError(f"Failed to copy multi_config.csv: {e}")

        # Verify the file was copied
        if not config_dest.exists():
            raise FileNotFoundError(f"multi_config.csv was not copied to {config_dest}")

        # Debug: Print the copied file before patching
        print("Copied multi_config.csv before patching:")
        with open(config_dest, "r") as file:
            print(file.read())

        # Patch the multi_config.csv file
        print("Patching multi_config.csv...")
        with open(config_dest, "r") as file:
            lines = file.readlines()

        patched_lines = []
        in_gene_expression = False
        section_replaced = False

        for line in lines:
            stripped = line.strip()

            if stripped == "[gene-expression]":
                # Start of gene-expression section
                in_gene_expression = True
                continue  # Skip this header for now

            # If another section begins, insert new gene-expression section before continuing
            if stripped.startswith("[") and in_gene_expression and not section_replaced:
                patched_lines.append("[gene-expression]\n")
                patched_lines.append("create-bam,true\n")
                patched_lines.append(f"reference,{self.ref_genome}\n")
                patched_lines.append(f"probe-set,{self.probe_file}\n")
                section_replaced = True
                in_gene_expression = False

            if not in_gene_expression:
                patched_lines.append(line)

        # Edge case: if gene-expression was last section and never replaced
        if in_gene_expression and not section_replaced:
            patched_lines.append("[gene-expression]\n")
            patched_lines.append("create-bam,true\n")
            patched_lines.append(f"reference,{self.ref_genome}\n")
            patched_lines.append(f"probe-set,{self.probe_file}\n")

        with open(config_dest, "w") as file:
            file.writelines(patched_lines)

        # Debug: Print the patched file
        print("Patched multi_config.csv:")
        with open(config_dest, "r") as file:
            print(file.read())




    def run_cellranger_multi(self):
        """Run the Cell Ranger multi command with module loading."""
        print("🚀 Running Cell Ranger multi with module loading...")

        # Define the module commands
        module_commands = (
            "set +u && "
            "module purge && "
            "module load Bioinformatics cellranger && "
            "module load snakemake && "
            "set -u"
        )

        # Define the Cell Ranger command
        output_dir = self.results_dir
        config_file = self.results_dir / "multi_config.csv"  # Use the patched version

        # Ensure the logs directory exists
        self.logs_dir.mkdir(parents=True, exist_ok=True)

        # Delete the output directory if it exists
        if output_dir.exists():
            print(f"⚠️ Output directory {output_dir} already exists. Deleting it to avoid conflicts.")
            shutil.rmtree(output_dir)

        # Verify the multi_config.csv file exists
        if not Path(config_file).exists():
            raise FileNotFoundError(f"multi_config.csv not found at {config_file}")

        # Properly quote the path to handle spaces
        quoted_config_file = f'"{config_file}"'
        quoted_output_dir = f'"{output_dir}"'

        cellranger_command = (
            f"cellranger multi --id={self.output_id} "
            f"--csv={quoted_config_file} "
            f"--output-dir={quoted_output_dir}"
        )

        # Combine the commands
        full_command = f"{module_commands} && {cellranger_command}"

        print(f"Executing command: {full_command}")
        subprocess.run(full_command, shell=True, check=True)

    def run(self):
        """Execute the Cell Ranger workflow."""
        print(f"🔬 Starting {self.workflow_name}...")
        self.validate_paths()
        self.prepare_multi_config()
        self.run_cellranger_multi()
        print(f"✅ {self.workflow_name} completed successfully!")