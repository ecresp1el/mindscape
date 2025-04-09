###############################################################################
# Cell Ranger Workflow using Snakemake
#
# This workflow runs `cellranger multi` using a config file copied into
# the test folder (/nfs/turbo/umms-parent/Manny_test).
#
# This avoids modifying the local project and ensures reproducible outputs.
#
# HOW TO RUN THIS WORKFLOW:
#   bash scripts/run_cellranger_pipeline.sh
###############################################################################

# Define expected final output to trigger the workflow
rule all:
    input:
        "10496-MW-reanalysis"

# Rule to run cellranger multi from a CSV config file
rule cellranger_multi:
    input:
        csv="/nfs/turbo/umms-parent/Manny_test/multi_config.csv"
    output:
        directory("10496-MW-reanalysis")
    shell:
        """
        cellranger multi \\
            --id=10496-MW-reanalysis \\
            --csv={input.csv}
        """
