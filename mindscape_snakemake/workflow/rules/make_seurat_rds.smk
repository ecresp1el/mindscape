import os

def get_matrix_dir(sample):
    return os.path.join(
        config["cellranger_base_dir"],
        sample,
        "count",
        "sample_filtered_feature_bc_matrix"
    )

def get_rds_output_path(sample):
    return os.path.join(
        config["project_output_dir"],
        "seurat_rds",
        f"{sample}.rds"
    )

rule init_output_dir:
    output:
        "init.done"
    run:
        os.makedirs(config["project_output_dir"], exist_ok=True)
        with open(output[0], "w") as f:
            f.write("initialized\n")

# Register a Seurat RDS rule for each sample
for sample in config["samples"]:
    rule_name = f"make_seurat_rds_{sample.replace('-', '_')}"
    rule:
        name: rule_name
        input:
            matrix_dir = get_matrix_dir(sample),
            init_flag = rules.init_output_dir.output
        output:
            rds = get_rds_output_path(sample)
        params:
            sample = sample
        script:
            "../scripts/create_seurat_rds.R"

# Final rule to trigger everything
rule all:
    input:
        [get_rds_output_path(sample) for sample in config["samples"]]