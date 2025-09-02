###############################################################################
# Cell Ranger Workflow (self-contained)
#
# This Snakefile is placed inside scripts/cellranger and is referenced by
# scripts/cellranger/scripts/create_project_cellranger_div90.sh.
#
# It assumes the wrapper has copied a patched multi_config.csv into the
# Snakemake working directory (the TEST_DIR), and that the Snakemake
# invocation specifies the desired Cell Ranger --id as the target.
#
# Usage (normally called via wrapper):
#   snakemake --directory "$TEST_DIR" --snakefile this_file.smk "div90-reanalysis"
#
# Notes:
# - The single target is a directory named after the Cell Ranger --id.
# - The wrapper prepares multi_config.csv in TEST_DIR and invokes this rule.
###############################################################################

rule cellranger_multi:
    input:
        # The wrapper prepares this CSV in the working directory
        csv="multi_config.csv"
    output:
        # The target directory is the Cell Ranger --id (wildcards.id)
        directory("{id}")
    shell:
        """
        cellranger multi \
            --id {wildcards.id} \
            --csv {input.csv}
        """
