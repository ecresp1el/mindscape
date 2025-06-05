rule sync_data:
    input:
        marker="results/create_project.done",
        config="config/config.yaml"
    output:
        marker="results/sync_data.done"
    script:
        "workflow/scripts/sync_ventral_div90_data.py"