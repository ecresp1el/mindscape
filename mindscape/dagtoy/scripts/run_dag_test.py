from mindscape.dagtoy.create_project.new import create_new_project

# Replace hardcoded path with dynamic creation
project_path, config_path = create_new_project(
    base_dir="mindscape/dagtoy", 
    project_name="TestDAGRun", 
    experimenter="Manny"
)

run_mindscape_flow(config_path)