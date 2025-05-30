from .base_workflow import BaseWorkflow
from pathlib import Path
import shutil
import subprocess

class CellRangerWorkflow(BaseWorkflow):
    def __init__(self, config):
        super().__init__(config)
        self.workflow_name = "CellRangerWorkflow"

        # Hardcoded parameters for now
        self.test_dir = Path(self.config["project_path"]) / "data"
        self.multi_config_source = "/nfs/turbo/umms-parent/Accessible_multi-config_csvs/multi-config.csv"
        self.probe_file = "/nfs/turbo/umms-parent/10X_Human_Refs/Probe_Set.csv"
        self.ref_genome = "/nfs/turbo/umms-parent/10X_Human_Refs/Refdata/refdata-gex-GRCh38-2020-A"
        self.output_id = "10496-MW-reanalysis"

    def prepare_multi_config(self):
        """Prepare the multi_config.csv file."""
        config_dest = self.test_dir / "multi_config.csv"

        # Ensure the reference genome exists
        if not Path(self.ref_genome).exists():
            raise FileNotFoundError(f"Reference genome not found at {self.ref_genome}")
        print(f"🧬 Using reference genome: {self.ref_genome}")

        # Copy the multi_config.csv file
        print(f"📄 Copying multi_config.csv to {config_dest}")
        shutil.copy(self.multi_config_source, config_dest)

        # Patch the multi_config.csv file
        print("🛠 Patching multi_config.csv...")
        with open(config_dest, "r") as file:
            lines = file.readlines()

        with open(config_dest, "w") as file:
            for line in lines:
                file.write(line)
                if "[gene-expression]" in line:
                    file.write("create-bam,true\n")
            file.write(f"reference,{self.ref_genome}\n")
            file.write(f"probe-set,{self.probe_file}\n")

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
        cellranger_command = (
            f"cellranger multi --id={self.output_id} "
            f"--csv={self.test_dir / 'multi_config.csv'}"
        )

        # Combine the commands
        full_command = f"{module_commands} && {cellranger_command}"

        print(f"Executing command: {full_command}")
        subprocess.run(full_command, shell=True, check=True)

    def run(self):
        """Execute the Cell Ranger workflow."""
        print(f"🔬 Starting {self.workflow_name}...")
        self.prepare_multi_config()
        self.run_cellranger_multi()
        print(f"✅ {self.workflow_name} completed successfully!")