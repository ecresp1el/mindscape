rule sync_ventral_div90_data:
    input:
        config = "config/config.yaml",
        trigger = "results/create_project.done"
    output:
        "results/sync_ventral_data.done"
    script:
        "workflow/scripts/sync_ventral_div90_data.py"