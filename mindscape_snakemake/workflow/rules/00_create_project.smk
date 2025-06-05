rule create_project:
    input:
        config="config/config.yaml"
    output:
        marker="results/create_project.done"
    script:
        "workflow/scripts/create_project.py"