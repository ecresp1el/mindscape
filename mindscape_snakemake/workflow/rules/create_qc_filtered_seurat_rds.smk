import os

def get_qc_filtered_rds_path(sample):
    return os.path.join(
        config["project_output_dir"],
        "qc_filtered",
        f"{sample}_qc_filtered.rds"
    )

def get_ridgeplot_path(sample):
    return os.path.join(
        config["project_output_dir"],
        "RidgePlots_QC_nFeature_nCounts_percent_Mito_sample",
        f"{sample}_QC_ridgeplot.png"
    )

for sample in config["samples"]:
    rule_name = f"create_qc_filtered_seurat_rds_{sample.replace('-', '_')}"
    rule:
        name: rule_name
        input:
            raw_rds = get_rds_output_path(sample)
        output:
            qc_rds = get_qc_filtered_rds_path(sample),
            ridgeplot = get_ridgeplot_path(sample)
        params:
            sample = sample
        resources:
            time = "1:00:00",
            mem_mb = 1500000,
            cpus = 16
        script:
            "../scripts/create_qc_filtered_seurat_rds.R"