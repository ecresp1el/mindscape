###############################################################################
# Cell Ranger Workflow using Snakemake
#
# This workflow runs `cellranger count` on one or more samples, using
# data and output paths defined in config/config.yaml.
#
# Designed specifically for the University of Michigan Great Lakes HPC.
#
# USERS SHOULD:
# - Only modify config/config.yaml to add samples or change paths
# - Not modify this file unless changing logic
#
# HOW TO RUN THIS:
#   snakemake -j 1 --snakefile mindscape/workflows/cellranger.smk --configfile config/config.yaml
###############################################################################

# Point to external configuration file (user-editable)
configfile: "config/config.yaml"

# Define expected final outputs for each sample
rule all:
    input:
        expand("{output_dir}/{sample}/outs/web_summary.html",
               output_dir=config["output_dir"],
               sample=config["samples"].keys())

# Rule to run Cell Ranger on a single sample
rule cellranger_count:
    output:
        "{output_dir}/{sample}/outs/web_summary.html"

    params:
        sample=lambda wc: wc.sample,
        fastqs=lambda wc: config["samples"][wc.sample]["fastqs"],
        transcriptome=lambda wc: config["samples"][wc.sample]["transcriptome"],
        outdir=lambda wc: os.path.join(config["output_dir"], wc.sample)

    shell:
        """
        module load python3.11-anaconda/2024.02
        module load Bioinformatics
        module load cellranger/9.0.1
        module load snakemake/7.32.4

        cellranger count \\
            --id={params.sample} \\
            --fastqs={params.fastqs} \\
            --sample={params.sample} \\
            --transcriptome={params.transcriptome} \\
            --localcores=8 \\
            --localmem=64 \\
            --output-dir={params.outdir}
        """
