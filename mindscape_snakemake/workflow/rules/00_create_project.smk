rule create_project:
    input:
        config="config/config.yaml"
    output:
        marker="results/create_project.done"
    script:
        "../scripts/create_project.py"